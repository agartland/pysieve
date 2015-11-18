"""
distance.py
Distance functions to be called by the analysis functions.

TODO:
 - Standardize inputs to distance functions to make the architecture more "plug-n-play"
 - Release a version of seqtools with the required functions to satisfy dependencies

Generally a distance function should have the following inputs:
 - Sieve data object
 - Parameters
   (Q: I use a dict of params to simplify the arguments, but is it confusing that each distance function requires different params?)

and should return either:
 - A 2D [nPTIDs x nSites] pd.DataFrame
 - A 3D [nPTIDs x nSites x nParams] numpy.ndarray for multiple parameter values


NOTES:
Currently it looks as though with low nPerms and high simulation repeats prepareBA and binding distance computations
are the speed bottlenck. I could probably speed up the distance comp by using the 4D MATLAB style BA matrix and numexpr
since this would vectorize the <, >, and sum operations."""

from seqdistance.matrices import binarySubst, addGapScores
from seqdistance import hamming_distance
import numpy as np
import pandas as pd
import itertools
from HLAPredCache import *
from decimal import Decimal

__all__=['_vxmatch_distance',
         '_binding_escape_distance',
         '_nepitope_distance',
         '_indel_escape_distance',
         '_prepareSeqBA',
         '_prepareAlignBA',
         '_prepareBTBA',
         '_prepareBA',
         '_similarity_score',
         '_findall',
         '_binding_scan_distance',
         '_indel_scan_distance',
         '_epitope_mismatch_distance',
         '_relative_binding_escape_distance']

def _vxmatch_distance(insertSeq, seqDf, params):
    """Computes a vaccine insert distance for each breakthrough sequence
    based on the specified "subst" similarity matrix.

    Used by the global and local vxmatch sieve analyses.

    NOTE:    
    Normalization is performed using a mean similarity across the sequence, but applied site-wise.
    This is important for stability in the normalization, but may be under estimating the contribution of
    rare amino-acids when non-binary and non-normalized similarity matrices are used (e.g. HIV PAM)

    Parameters
    ----------
    insertSeq : str
        Amino acid sequence of the vaccine insert/immunogen
    seqDf : pd.DataFrame
        Contains the breakthrough sequences (see pysieve.data for details)
    params : dict
        Should contain a similarity matrix 'subst', a dict of amino acid similarity scores.

    Returns
    -------
    dist : pd.DataFrame
        Distance matrix [ptids x sites] with PTID row-index and sites (0-indexed) as columns."""

    try:
        subst = params['subst']
    except KeyError:
        print 'Using default binary substitution matrix.'
        subst = addGapScores(binarySubst, None)

    N = seqDf.shape[0]
    nSites = len(insertSeq)

    sim = np.ones((N, nSites))
    
    """Compute insert2insert similarity for normalization.
    Similarity of insert to itself is high (exactly 1 for a hamming distance)"""
    insert2insert = np.nanmean(_similarity_score(insertSeq, insertSeq, subst = subst))
    
    for ptidi,ptid in enumerate(seqDf.index):
        """Similarity of a sequence to itself is also high and will be used in the normalization"""
        seq2seq = np.nanmean(_similarity_score(seqDf.seq[ptid], seqDf.seq[ptid], subst = subst))
        sim[ptidi,:] = _similarity_score(seqDf.seq[ptid], insertSeq, subst = subst, denominator = seq2seq + insert2insert)
    
    """similarity --> distance only works with correct normalization"""
    dist = 1 - sim

    return pd.DataFrame(dist, index = seqDf.index, columns = np.arange(nSites))

def _nepitope_distance(insertBA, btBA, params):
    """Creates a distance matrix (DataFrame) [N x sites]
    indicating the number of PTIDs predicted to have an insert epitope for each kmer

    Parameters
    ----------
    insertBA : pd.DataFrame
        Row index as HLAs and colums as sites, shape [len(uHLA4) x nSites]
    btBA : dict
        Dict contains keys for (1) "validHLAs" HLAs used, (2) "ba" array([nHLAs (4 typically), nSites]) and (3) "ptid"
    params : dict
        Should contain 'binding' and 'nmer' parameters.

    Returns
    -------
    dist : pd.DataFrame
        Distance matrix [ptids x sites] with PTID row-index and sites (kmer start positions 0-indexed) as columns."""

    N = len(btBA['ptid'])
    nSites = insertBA.shape[1]
    
    dist = np.nan * np.ones((N,nSites))
    for ptidi,ptid,ba,h in zip(np.arange(N), btBA['ptid'], btBA['ba'], btBA['validHLAs']):
        tmpInsert = insertBA.loc[h]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(h, return_index = True)

        """Note that the original binding escape distance sums across HLAs not any's """
        #dist[ptidi,:]=np.squeeze(np.any((tmpInsert<params['binding']),axis=0))
        dist[ptidi,:] = np.squeeze(np.sum((tmpInsert < params['binding']).values[uniqi,:], axis = 0))

    return pd.DataFrame(dist, index = btBA['ptid'], columns = np.arange(nSites))

def _indel_escape_distance(insertSeq, insertBA, seqDf, btBA, params):
    """Creates a distance matrix (DataFrame) [N x sites] with PTID rows and sites as columns
    populated with the HLA binding escape count distance

    Parameters
    ----------
    insertSeq : str
        Amino acid sequence of the vaccine insert/immunogen
    insertBA : pd.DataFrame
        Row index as HLAs and colums as sites, shape [len(uHLA4) x nSites]
    seqDf : pd.DataFrame
        Contains the breakthrough sequences (see pysieve.data for details)
    btBA : dict
        Dict contains keys for (1) "validHLAs" HLAs used, (2) "ba" array([nHLAs (4 typically), nSites]) and (3) "ptid"
    params : dict
        Should contain 'binding' and 'nmer' parameters.

    Returns
    -------
    dist : pd.DataFrame
        Distance matrix [ptids x sites] with PTID row-index and sites (kmer start positions 0-indexed) as columns."""
    
    N = seqDf.shape[0]
    nSites = insertBA.shape[1]
    
    dist = np.nan * np.ones((N, nSites))
    for ptidi,ptid in enumerate(seqDf.index):
        """An indel is 'shared' if its found in both the insert and breakthrough kmer"""
        unsharedIndel = np.zeros((len(btBA['validHLAs'][ptidi]), nSites), dtype = bool)
        for sitei in range(nSites):
            """grabKmer returns None if the kmer is invalid (e.g. off the end of the sequence)"""
            insertMer,_nonGapped = grabKmer(insertSeq, sitei, params['nmer'])
            btMer,_nonGapped = grabKmer(seqDf.seq[ptid], sitei, params['nmer'])
            
            """If the insert and the bt mer don't have gaps in the same places then there is an indel!"""
            if not insertMer is None:
                if btMer is None or not np.all(np.array(_findall(insertMer,'-')) == np.array(_findall(btMer,'-'))):
                    unsharedIndel[:,sitei] = True
      
        tmpInsert = insertBA.loc[btBA['validHLAs'][ptidi]]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(btBA['validHLAs'][ptidi], return_index=True)

        """ANY HLA: Escape count is 1 if the kmer binds and there is an indel in the BT seq"""
        #dist[ptidi,:] = np.squeeze(np.any((tmpInsert < params['binding']) & unsharedIndel, axis=0))

        """SUM HLAs: count 1 escape per HLA that binds (not including homozygous alleles"""
        dist[ptidi,:] = np.squeeze(np.sum(((tmpInsert < params['binding']).values & unsharedIndel)[uniqi,:],axis=0))
        
    return pd.DataFrame(dist, index = btBA['ptid'], columns = np.arange(nSites))

def _binding_escape_distance(insertSeq, insertBA, seqDf, btBA, params):
    """Creates a distance matrix (pd.DataFrame) [ptids x sites]
    populated with the HLA binding escape count distance.

    TODO:
     - Handle 15mer 'unique core' distances
     - Standardize the input arguments to match other distances
       (or at least other T cell based distances)
     - Make the handling of homozygous alleles and multiple escapes per kmer
       a parameter.

    Parameters
    ----------
    insertSeq : str
        Amino acid sequence of the vaccine insert/immunogen
    insertBA : pd.DataFrame
        Row index as HLAs and colums as sites, shape [len(uHLA4) x nSites]
    seqDf : pd.DataFrame
        Contains the breakthrough sequences (see pysieve.data for details)
    btBA : dict
        Dict contains keys for (1) "validHLAs" HLAs used, (2) "ba" array([nHLAs (4 typically), nSites]) and (3) "ptid"
    params : dict
        Should contain 'binding', 'escape' and 'nmer' parameters.

    Returns
    -------
    dist : pd.DataFrame
        Distance matrix [ptids x sites] with PTID row-index and sites (kmer start positions 0-indexed) as columns."""

    N = len(btBA['ptid'])
    nSites = insertBA.shape[1]
    
    """Don't count a binding escape if there's also an insertion/deletion there
    (these distances should be mutually exclusive)
    Import to make indelDist 0s and 1s to work for this purpose"""
    
    indelDist = (_indel_escape_distance(insertSeq, insertBA, seqDf, btBA, params).values > 0).astype(int)

    dist = np.nan * np.ones((N,nSites))
    for ptidi,ptid,ba,h in zip(np.arange(N), btBA['ptid'], btBA['ba'], btBA['validHLAs']):
        """Maxtrix of binding affinities for the hla alleles in h"""
        tmpInsert = insertBA.loc[h]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(h, return_index=True)

        """ANY HLA: For each HLA (typically 4 per PTID), if it meets the criteria for this kmer then its an escape"""
        #dist[ptidi,:] = np.squeeze(np.any((tmpInsert<params['binding']) & (ba>params['escape']),axis=0)) * (1-indelDist[ptidi,:])

        """SUM HLAS: Count multiple escapes per kmer if the person has multiple alleles with escape"""
        dist[ptidi,:] = np.squeeze(np.sum(((tmpInsert < params['binding']).values & (ba > params['escape']))[uniqi,:], axis=0)) * (1 - indelDist[ptidi,:])

    return pd.DataFrame(dist, index=btBA['ptid'], columns=np.arange(nSites))

def _prepareSeqBA(seq, hlas, ba, k, ignoreGappedKmers=False, getFast=False):
    """Prepare a matrix of binding affinities for all kmers in seq and all hlas.

    Parameters
    ----------
    seqs : collection
        Aligned sequences/strings
    hlas : collection
        HLA alleles
    ba : dict/HLAPredCache
        Contains all neccessary predicted binding affinities for propagating the matrix.
        Keys are tuples (allele, peptide)
    k : int
        Length of the peptides.
    ignoreGappedKmers : bool
        If False then kmer continues for full k AA,
        but if True then throws out all kmers with a gap at any location in any bt sequence.
    getFast : bool
        If True, uses the getFast method of hlaPredCache, w/o error checking

    Returns
    -------
    baMat : pd.DataFrame, [nHLAs x nSites]
        Matrix of binding affinities with rows HLAs and columns as 0-indexed kmer start positions."""
    
    nSites = len(seq)
    baMat = np.nan * np.ones((len(hlas),nSites))

    if getFast:
        baFunc = lambda t: ba.getFast(t)
        """Replace all '*' here just in case, if trying to use getFast"""
        originalHLAs = hlas
        hlas = [h.replace('*','_') for h in hlas]
    else:
        originalHLAs = hlas
        baFunc = lambda t: ba[t]

    for sitei in range(nSites):
        """ngmer is None if the insert kmer starts with a gap '-', leave these as nan"""
        gmer,ngmer = grabKmer(seq,sitei,k)
        if not ignoreGappedKmers:
            """Use ngmer which starts at sitei and grabs the next nmer AAs (not counting gaps)"""
            mer = ngmer
        else:
            mer = gmer

        if (not mer is None) and (not '-' in mer):
            for hlai,hla in enumerate(hlas):
                if isvalidHLA(hla):
                    baMat[hlai,sitei] = baFunc((hla,mer))

    return pd.DataFrame(baMat, index = originalHLAs, columns = np.arange(nSites))

def _prepareAlignBA(seqs, hlas, ba, k, ignoreGappedKmers = False, getFast = False):
    """Prepare a matrix of binding affinities for all kmers, all HLAs and all seqs.

    Parameters
    ----------
    seqs : collection
        Aligned sequences/strings
    hlas : collection
        HLA alleles
    ba : dict/HLAPredCache
        Contains all neccessary predicted binding affinities for propagating the matrix.
        Keys are tuples (allele, peptide)
    k : int
        Length of the peptides.
    ignoreGappedKmers : bool
        If False then kmer continues for full k AA,
        but if True then throws out all kmers with a gap at any location in any bt sequence.
    getFast : bool
        If True, uses the getFast method of hlaPredCache, w/o error checking

    Returns
    -------
    baMat : ndarray [nSeqs x nHLAs x nSites]
        Matrix of binding affinities"""
    
    nSites = int(np.median([len(s) for s in seqs]))
    baMat = np.nan * np.ones((len(seqs), len(hlas), nSites))

    """Go through each person, creating a bt BA [len(hla) x nSites] and assign to the big matrix"""
    for seqi,seq in enumerate(seqs):
        baMat[seqi,:,:] = _prepareSeqBA(seq, hlas, ba, k, ignoreGappedKmers, getFast).values

    """Ignore any kmer that had a gap in any BT sequence"""
    if ignoreGappedKmers:
        """If the BA is nan for all HLAs then the kmer must have had a gap.
        If aross people, any kmer had a gap then set all people nan there"""
        badSites = np.any(np.all(np.isnan(baMat), axis=1), axis=0)
        """TODO: this could be simplified using numpy broadcasting"""
        baMat[np.tile(badSites[None,None,:], (baMat.shape[0], baMat.shape[1],1))] = np.nan

    return baMat

def _prepareBTBA(data, ba, params):
    """Prepare matrix of log-IC50 binding affinities of BREAKTHROUGH sequences given the dict/HLAPredCache ba

    Only the BA for alleles expressed by each subject are returned.

    Q:What calls this function and why?

    Parameters
    ----------
    data : pysieve.sieveData object
    ba : dict/hlaPredCache
    params : dict
        Required keys: nmer, ignoreGappedKmers, getFast
    
    Returns
    -------
    btBA : dict of lists with rows per ptid
        Dict keys for (1) "validHLAs" HLAs used,
                      (2) "ba" array([nHLAs (4 typically), nSites])
                      (3) "ptid"
    """
    nSites = len(data.insertSeq)

    fullBTBA = _prepareAlignBA(data.seqDf.seq,
                               data.uHLA4,
                               ba,
                               params['nmer'],
                               ignoreGappedKmers = params['ignoreGappedKmers'],
                               getFast = params['getFast'])

    """Go through each person, creating a bt BA [len(hla) x nSites]"""
    btBA = {'ptid':[],'validHLAs':[],'ba':[]}
    for ptid,row in data.seqDf.iterrows():
        ptidi = list(data.seqDf.index).index(ptid)
        btSeq = row['seq']
        HLAs = data.ptidDf.hla[ptid]
        validHLAs = []
        for hla in HLAs:
            if isvalidHLA(hla):
                validHLAs.append(hla)
        
        """New way using the full BTBA"""
        tmpba = np.nan * np.ones((len(validHLAs),nSites))
        for i,h in enumerate(validHLAs):
            hlai = list(data.uHLA4).index(h)
            tmpba[i,:]=fullBTBA[ptidi,hlai,:]
        
        btBA['ptid'].append(ptid)
        btBA['validHLAs'].append(validHLAs)
        btBA['ba'].append(tmpba)
    
    return btBA

def _prepareBA(data, ba, params):
    """Prepare matrices of log-IC50 HLA binding affinities for insert and breakthrough sequences
    Used by several T-cell based sieve distances.

    TODO:
    Currently this function returns two different btBA data objects depending on fullBT
    This is probably not a good idea.
    
    Parameters
    ----------
    data : pysieve.sieveData object
    ba : dict/hlaPredCache
    params : dict
        Required keys: nmer, fullBT, ignoreGappedKmers, getFast
    
    Returns
    -------
    insertBA : pd.DataFrame, shape [len(uHLA4) x nSites]
        Row index as HLAs and colums as start positions
            
    btBA : variable
        Dict of lists with rows per ptid, dict keys for
            (1) "validHLAs" HLAs used,
            (2) "ba" array([nHLAs (4 typically), nSites])
            (3) "ptid"
        OR if fullBT is True, a 3D ndarray [nSeqs x nHLAs x nSites]
    method : str
        Describing the prediction method used (from ba)"""

    fullBT = params.get('fullBT',False)
    getFast = params.get('getFast',False)
    ignoreGappedKmers = params.get('ignoreGappedKmers',False)

    """Create matrix of insert BA [len(uHLA4) x nSites]"""
    insertBA = _prepareSeqBA(data.insertSeq,
                             data.uHLA4,
                             ba,
                             params['nmer'],
                             ignoreGappedKmers = params['ignoreGappedKmers'],
                             getFast = params['getFast'])

    if not fullBT:
        btBA = _prepareBTBA(data,ba,params)
    else:
        """New structure required by _relative_binding_escape distance"""
        btBA = _prepareAlignBA(data.seqDf.seq,
                               data.uHLA4,
                               ba,
                               params['nmer'],
                               ignoreGappedKmers = params['ignoreGappedKmers'],
                               getFast = params['getFast'])
    
    return insertBA, btBA, ba.predictionMethod

def _similarity_score(seq1, seq2, subst, denominator = 2):
    """Return a vector of site-wise similarities between two sequences based on a substitution matrix (dict).
    
    Optionally can give a denominator for normalization.
    Example denominator: sim11 + sim22 which is the sum of seq1 to itself and seq2 to itself.
    
    Denominator can be supplied as a vector, in which case the normalization is done site-wise or
    as a scalar in which case it is equivalent to applying the normalization for the whole sequence
    (even though this function will still give the distance site-wise)
    
    By default there is no normalization (denominator = 2). 
    This can create a problem for vxmatch using similarity matrices other than binarySubst"""

    """Similarity between seq1 and seq2 using the substitution matrix subst"""
    sim12 = np.array([i for i in itertools.imap(lambda a,b: subst.get((a,b), subst.get((b,a))), seq1, seq2)])

    return (2 * sim12) / denominator

def _findall(s, item):
    """Return index of each element equal to item in s"""
    return [i for i,x in enumerate(s) if x == item]

def _indel_scan_distance(insertSeq, insertBA, seqDf, btBA, nmer, paramPairs, minDelta):
    """Creates a distance matrix ndarray [N x sites x params]
    populated with the indel escape count distance"""

    N = seqDf.shape[0]
    nSites = insertBA.shape[1]
    nParams = len(paramPairs)
    minDelta = Decimal('%1.1f' % minDelta)
    
    dist = np.nan*np.ones((N,nSites,nParams))
    for ptidi,ptid in enumerate(seqDf.index):
        unsharedIndel = np.zeros((len(btBA['validHLAs'][ptidi]),nSites), dtype=bool)
        """Determine if there are any 'unshared' gaps in each kmer"""
        for sitei in xrange(nSites):
            insertMer,_nonGapped = grabKmer(insertSeq,sitei,nmer)
            btMer,_nonGapped = grabKmer(seqDf.seq[ptid],sitei,nmer)
            """If the insert and the bt mer don't have gaps in the same places then there is an indel!"""
            if not insertMer is None:
                if btMer is None or not np.all(np.array(_findall(insertMer,'-')) == np.array(_findall(btMer,'-'))):
                    unsharedIndel[:,sitei] = True

        tmpInsert = insertBA.loc[btBA['validHLAs'][ptidi]]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(btBA['validHLAs'][ptidi], return_index=True)

        """For each parameter pair compute the distance"""
        for parami,pp in enumerate(paramPairs):
            """A pair of binding and escape thresholds are only valid if escape - binding > delta"""
            if (pp[1]-pp[0]) > minDelta:
                """If any of the HLAs (typically 4 per PTID) meets the criteria for this kmer then its an escape"""
                #dist[ptidi,:,parami]=np.squeeze(np.any((tmpInsert < pp[0]) & unsharedIndel,axis=0))
                """The original binding escape distance, for EACH of the HLAs that meets the criteria...its an escape"""
                dist[ptidi,:,parami] = np.squeeze(np.sum(((tmpInsert<np.float64(pp[0])).values & unsharedIndel)[uniqi,:],axis=0))
    return dist

def _binding_scan_distance(insertSeq, insertBA, seqDf, btBA, paramPairs, minDelta, nmer):
    """Creates a distance matrix ndarray [N x sites x params]
    populated with the HLA binding escape count distance"""
    N = len(btBA['ptid'])
    nSites = insertBA.shape[1]
    nParams = len(paramPairs)
    minDelta = Decimal('%1.1f' % minDelta)
    
    """Don't count a binding escape if there's also an indel there (these distances should be mutually exclusive)
    Import to dichotomize this distance when used here though"""
    #indelDist = (_indel_scan_distance(insertSeq, insertBA, seqDf, btBA, nmer, paramPairs, minDelta) > 0).astype(np.int64)

    """Nan is default value!
    This means that the following will be nan in dist:
        non-existent insert kmers (ie start with '-') (ACTUALLY THESE SHOULD BE ZEROS)
        invalid parameters (insufficient delta)
        filteredDist will have nans at kmers that have been filtered out
        it doesn't appear that there are other nans
    This is important because it means in the compstats for binding escape i could set nan to 0
    """
    dist = np.nan * np.ones((N, nSites, nParams))
    dummydist = np.nan * np.ones((N, nSites, nParams))
    for ptidi,ptid,ba,h in zip(np.arange(N), btBA['ptid'], btBA['ba'], btBA['validHLAs']):
        tmpInsert = insertBA.loc[h]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(h, return_index=True)

        for parami,pp in enumerate(paramPairs):
            """A pair of binding and escape thresholds are only valid if escape - binding > delta"""
            if (pp[1]-pp[0]) > minDelta:
                tmp=(tmpInsert < np.float64(pp[0])).values & (ba > np.float64(pp[1]))
                """Sum across HLAs"""
                dist[ptidi,:,parami] = np.sum(tmp[uniqi,:], axis=0) #* (1-indelDist[ptidi,:,parami])
                dummydist[ptidi,:,parami] = 1
    return dist

def _epitope_mismatch_distance(seqDf, insertSeq, insertDf, insertBA, btBA, params):
    """Creates a distance matrix (DataFrame) [N x sites]
    indicating the PTIDs with an insert epitope and a
    breakthrough with greater than mmTolerance substitutions relative to the reference.

    Parameters
    ----------
    seqDf : pd.DataFrame
        PTID as index with seq column containing breakthrough sequences
    insertSeq : str
        Amino acid sequence of the insert
    insertBA : pd.DataFrame
        Row index as HLAs and colums as sites, shape [len(uHLA4) x nSites]
    btBA : dict
        Dict contains keys for (1) "validHLAs" HLAs used, (2) "ba" array([nHLAs (4 typically), nSites]) and (3) "ptid"
    params : dict
        Should contain binding, nmer, ignoreGappedKmers, and mmTolerance parameters.

    Returns
    -------
    dist : pd.DataFrame
        Distance matrix [ptids x sites] with PTID row-index and sites (kmer start positions 0-indexed) as columns."""

    N = len(btBA['ptid'])
    nSites = insertBA.shape[1]

    mmCount = np.zeros((N, nSites))
    for sitei in range(nSites):
        """ngmer is None if the insert kmer starts with a gap '-', leave these as nan"""
        igmer,ingmer = grabKmer(insertSeq, sitei, params['nmer'])
        if not params['ignoreGappedKmers']:
            """Use ngmer which starts at sitei and grabs the next nmer AAs (not counting gaps)"""
            imer = ingmer
        else:
            imer = igmer
        for ptidi,ptid in zip(np.arange(N), btBA['ptid']):
            btgmer,btngmer = grabKmer(seqDf.seq.loc[ptid], sitei, params['nmer'])
            if not params['ignoreGappedKmers']:
                """Use ngmer which starts at sitei and grabs the next nmer AAs (not counting gaps)"""
                btmer = btngmer
            else:
                btmer = btgmer
        mmCount[ptidi,sitei] = hamming_distance(imer, btmer, asStrings=True)
    
    dist = np.nan * np.ones((N, nSites))
    for ptidi,ptid,ba,h in zip(np.arange(N), btBA['ptid'], btBA['ba'], btBA['validHLAs']):
        tmpInsert = insertBA.loc[h]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = np.unique(h, return_index=True)

        """ANY HLA allele binds AND mmCount is greater than mmTolerance"""
        dist[ptidi,:] = np.squeeze(np.any((tmpInsert < params['binding']).values[uniqi,:], axis=0))
        dist[ptidi,:] = dist[ptidi,:] & (mmCount[ptidi,:] > params['mmTolerance'])

    return pd.DataFrame(dist.astype(np.float64), index=btBA['ptid'], columns=np.arange(nSites))

def _identifyNoPressureBT(insertBA, hlaMatrix, params):
    """Identify the bt kmers that aren't under pressure from potential vaccine epitopes (the person doesn't have any binding HLAs for the kmer location)
    Returns a numpy array [N x 1 x nSites]
    insertBA: [nHLAs x nSites] pd.DataFrame
    hlaMatrix: boolean matrix of HLA expression [ptid x nHLAs]
    params: binding"""
    N = hlaMatrix.shape[0]
    nSites = insertBA.shape[1]

    """Identify the binding HLAs at each site: [nHLAs x nSites] bool"""
    bindingHLAs = insertBA < params['binding']

    """Identify the breakthrough kmers (BA) associated with people lacking a binding allele at each site"""
    """[ptid x hla x site]"""
    tmpHLA = np.tile(hlaMatrix.values[:,:,None],(1,1,nSites))
    tmpBinding = np.tile(bindingHLAs.values[None,:,:],(N,1,1))

    """Use this to index into btBA to pull out the BA for kmers that were not under pressure by the HLA (e.g. btBA[:,hlai,:][noPressureBT])"""
    """[ptid x 1 x sites]"""
    """For each bt kmer, for it not to be under pressure,
    it must either not be an insert binder or not be an HLA allele that the person expresses, across all HLAs"""
    noPressureBT = np.all(~tmpBinding | ~tmpHLA, axis=1)[:,None,:]

    """If any sites don't have any bt kmers that are not under pressure then raise an exception"""
    if np.any(np.squeeze(noPressureBT.sum(axis=0)) == 0):
        raise Exception("Can't compute escape threshold for all kmers!")
    return noPressureBT


def _relative_binding_escape_distance(insertSeq,insertBA,seqDf,btBA,hlaMatrix,params,lookupDf):
    """Creates a distance matrix (DataFrame) [N x sites] with PTID rows and sites as columns
    populated with the HLA binding escape count distance
    
    This is "Allan's distance" and it differs from the binding_escape distance
    only in how the escape threshold is determined. The escape threshold is the
    median of the BA of the non-binding HLAs with a given BT peptide

    insertSeq: AA str
    insertBA: [nHLAs x nSites] pd.DataFrame
    seqDf: DataFrame with ptid rows and column seq containing BT seqs
    btBA: FULL btBA matrix [nSeqs x nHLAs x nSites] ndarray (neccessary for this method)
    hlaMatrix: boolean matrix of HLA expression [ptid x nHLAs]
    params: binding"""

    N = btBA.shape[0]
    nSites = insertBA.shape[1]

    """Don't count a binding escape if there's also an indel there (these distances should be mutually exclusive)
    Import to make indelDist 0s and 1s to work for this purpose"""
    #indelDist=(_indel_escape_distance(insertSeq,insertBA,seqDf,btBA,params).values > 0).astype(int64)

    """Identify the breakthrough kmers (BA) associated with people lacking a binding allele at each site"""
    noPressureBT = _identifyNoPressureBT(insertBA, hlaMatrix, params)
    
    """with open(DATA_PATH + 'STEP/andrew_rel_binding_escape.csv','w') as fh:
        print >> fh, 'position,seqid,ptid,hla,insert_peptide,bt_peptide,rbe,be,cutoff,ndiff,y,BET'"""

    dist = np.nan * np.ones((N,nSites))
    for ptidi,ptid in enumerate(seqDf.index):
        """Slice insertBA so that it only contains the valid HLAs of the current ptid"""
        validHLAs = [h for h in hlaMatrix.columns[hlaMatrix.ix[ptid]] if isvalidHLA(h)]
        validHLAInd = np.array([h in validHLAs for h in hlaMatrix.columns])
        """[4 x nSites]"""
        tmpInsert = insertBA.ix[validHLAs].values

        """Do not double count escapes from homozygous alleles"""
        uValidHLAs,uniqi = np.unique(validHLAs, return_index=True)

        """For each HLA (typically 4 per PTID), if it meets the criteria for this kmer then its an escape"""

        """[4 x nSites]"""
        insertBinders = (tmpInsert < params['binding'])

        """This is a complicated step:
            (1) Pull out the BA for the HLAs associated with this ptid [N x 4 x nSites]
            (2) Tile the noPressureBT index matrix along the hla axis [N x 4 nSites]
            (3) Index (1) using (2) to yield a matrix of the btBAs with the current HLAs (4) at the sequences of people that had no binding HLA
                [nonBinders x 4 x nSites]
            (4) Take the median across nonBinders
                [4 x nSites]
        """

        """Set all btBA that are under pressure from an HLA to nan prior to taking the median across bt kmers"""
        tmpNoPressureInd = np.tile(noPressureBT, (1,len(validHLAs),1))
        tmpBtBA = btBA[:,validHLAInd,:]
        tmpBtBA[~tmpNoPressureInd] = np.nan
        """tauThresholds [4 x nSites]"""
        tauThresholds = np.nanmedian(tmpBtBA, axis=0)

        """nonBTBinders [4 x nSites]"""
        nonBTBinders = (btBA[ptidi,validHLAInd,:] > tauThresholds) & (tmpInsert < btBA[ptidi,validHLAInd,:])
        escapes = (insertBinders & nonBTBinders)[uniqi,:]
        dist[ptidi,:] = np.squeeze(np.sum(escapes, axis=0))# * (1-indelDist[ptidi,:])

        """for hi in arange(escapes.shape[0]):
            hla=uValidHLAs[hi]
            hlai=insertBA.index==hla
            seqid=lookupDf.index[lookupDf.ptid==ptid][0]
            for kmeri in arange(escapes.shape[1]):
                gapped,insertKmer=grabKmer(insertSeq,kmeri,k=params['nmer'])
                gapped,btKmer=grabKmer(seqDf.seq[ptid],kmeri,k=params['nmer'])
                try:
                    hdist=hamming_distance(insertKmer,btKmer)
                except TypeError:
                    hdist=nan
                
                if seqid=='5021709' and kmeri==10 and hla=='A_2301':
                    continue
                    raise Exception()
                
                #print >> fh, 'position,seqid,ptid,hla,insert_peptide,bt_peptide,rbe,be,cutoff,ndiff,y,BET'
                print >> fh, '%d,%s,%s,%s,%s,%s,%1.2f,%1.2f,%1.2f,%1.0f,%d,%1.1f' % (kmeri+1,
                            seqid,ptid,hla,insertKmer,btKmer,insertBA[kmeri][hlai],btBA[ptidi,hlai,kmeri],
                            tauThresholds[uniqi,:][hi,kmeri],hdist,escapes[hi,kmeri],params['binding'])"""
    return pd.DataFrame(dist, index=seqDf.index, columns=insertBA.columns)