"""
distance.py
Distance functions to be called by the analysis functions, including:
    vxmatch_global
    vxmatch_site
    binding_escape_global
    binding_escape_kmer

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

from __future__ import division
from seqtools import *
from numpy import *
import pandas as pd
import itertools
from nanfunctions import *
from hla_prediction import *

__all__=['_vxmatch_distance',
         '_binding_escape_distance',
         '_indel_escape_distance',
         '_prepareSeqBA',
         '_prepareAlignBA',
         '_prepareBTBA',
         '_prepareBA',
         '_similarity_score']

def _vxmatch_distance(insertSeq,seqDf,params):
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

    sim = ones((N,nSites))
    
    """Compute insert2insert similarity for normalization.
    Similarity of insert to itself is high (exactly 1 for a hamming distance)"""
    insert2insert = nanmean(_similarity_score(insertSeq, insertSeq, subst = subst))
    
    for ptidi,ptid in enumerate(seqDf.index):
        """Similarity of a sequence to itself is also high and will be used in the normalization"""
        seq2seq = nanmean(_similarity_score(seqDf.seq[ptid], seqDf.seq[ptid], subst = subst))
        sim[ptidi,:] = _similarity_score(seqDf.seq[ptid], insertSeq, subst = subst, denominator = seq2seq + insert2insert)
    
    """similarity --> distance only works with correct normalization"""
    dist = 1 - sim

    return pd.DataFrame(dist, index = seqDf.index, columns = arange(nSites))

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
    
    dist = nan*ones((N,nSites))
    for ptidi,ptid in enumerate(seqDf.index):
        """An indel is 'shared' if its found in both the insert and breakthrough kmer"""
        unsharedIndel = zeros((len(btBA['validHLAs'][ptidi]), nSites), dtype = bool)
        for sitei in xrange(nSites):
            """grabKmer returns None if the kmer is invalid (e.g. off the end of the sequence)"""
            insertMer,_nonGapped = grabKmer(insertSeq, sitei, params['nmer'])
            btMer,_nonGapped = grabKmer(seqDf.seq[ptid], sitei, params['nmer'])
            
            """If the insert and the bt mer don't have gaps in the same places then there is an indel!"""
            if not insertMer is None:
                if btMer is None or not all(array(findall(insertMer,'-')) == array(findall(btMer,'-'))):
                    unsharedIndel[:,sitei] = True
      
        tmpInsert = insertBA.ix[btBA['validHLAs'][ptidi]]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = unique(btBA['validHLAs'][ptidi],return_index=True)

        """ANY HLA: Escape count is 1 if the kmer binds and there is an indel in the BT seq"""
        #dist[ptidi,:] = squeeze(any((tmpInsert < params['binding']) & unsharedIndel, axis=0))

        """SUM HLAs: count 1 escape per HLA that binds (not including homozygous alleles"""
        dist[ptidi,:] = squeeze(sum(((tmpInsert < params['binding']).values & unsharedIndel)[uniqi,:],axis=0))
        
    return pd.DataFrame(dist, index = btBA['ptid'], columns = arange(nSites))

def _binding_escape_distance(insertSeq,insertBA,seqDf,btBA,params):
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

    dist = nan * ones((N,nSites))
    for ptidi,ptid,ba,h in zip(arange(N), btBA['ptid'], btBA['ba'], btBA['validHLAs']):
        """Maxtrix of binding affinities for the hla alleles in h"""
        tmpInsert = insertBA.ix[h]

        """Do not double count escapes from homozygous alleles"""
        dummy,uniqi = unique(h, return_index = True)

        """ANY HLA: For each HLA (typically 4 per PTID), if it meets the criteria for this kmer then its an escape"""
        #dist[ptidi,:] = squeeze(any((tmpInsert<params['binding']) & (ba>params['escape']),axis=0)) * (1-indelDist[ptidi,:])

        """SUM HLAS: Count multiple escapes per kmer if the person has multiple alleles with escape"""
        dist[ptidi,:] = squeeze(sum(((tmpInsert < params['binding']).values & (ba > params['escape']))[uniqi,:], axis=0)) * (1 - indelDist[ptidi,:])

    return pd.DataFrame(dist, index = btBA['ptid'], columns = arange(nSites))

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
    baMat = nan * ones((len(hlas),nSites))

    if getFast:
        baFunc = lambda t: ba.getFast(t)
        """Replace all '*' here just in case, if trying to use getFast"""
        originalHLAs = hlas
        hlas = [h.replace('*','_') for h in hlas]
    else:
        originalHLAs = hlas
        baFunc = lambda t: ba[t]

    for sitei in xrange(nSites):
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

    return pd.DataFrame(baMat, index = originalHLAs, columns = arange(nSites))

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
    
    nSites = int(median([len(s) for s in seqs]))
    baMat = nan * ones((len(seqs), len(hlas), nSites))

    """Go through each person, creating a bt BA [len(hla) x nSites] and assign to the big matrix"""
    for seqi,seq in enumerate(seqs):
        baMat[seqi,:,:] = _prepareSeqBA(seq, hlas, ba, k, ignoreGappedKmers, getFast).values

    """Ignore any kmer that had a gap in any BT sequence"""
    if ignoreGappedKmers:
        """If the BA is nan for all HLAs then the kmer must have had a gap.
        If aross people, any kmer had a gap then set all people nan there"""
        badSites = any(all(isnan(baMat), axis=1), axis=0)
        """TODO: this could be simplified using numpy broadcasting"""
        baMat[tile(badSites[None,None,:], (baMat.shape[0], baMat.shape[1],1))] = nan

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

    fullBTBA = _prepareAlignBA(data.seqDf.seq, data.uHLA4, ba, params['nmer'], ignoreGappedKmers=params['ignoreGappedKmers'], getFast=params['getFast'])

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
        tmpba = nan*ones((len(validHLAs),nSites))
        for i,h in enumerate(validHLAs):
            hlai = list(data.uHLA4).index(h)
            tmpba[i,:]=fullBTBA[ptidi,hlai,:]
        
        btBA['ptid'].append(ptid)
        btBA['validHLAs'].append(validHLAs)
        btBA['ba'].append(tmpba)
    
    return btBA

def _prepareBA(data,ba,params):
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
    insertBA = _prepareSeqBA(data.insertSeq,data.uHLA4, ba, params['nmer'], ignoreGappedKmers=params['ignoreGappedKmers'], etFast=params['getFast'])

    if not fullBT:
        btBA = _prepareBTBA(data,ba,params)
    else:
        """New structure required by _relative_binding_escape distance"""
        btBA = _prepareAlignBA(data.seqDf.seq,data.uHLA4, ba, params['nmer'], ignoreGappedKmers=params['ignoreGappedKmers'], getFast=params['getFast'])
    
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
    sim12 = array([i for i in itertools.imap(lambda a,b: subst.get((a,b), subst.get((b,a))), seq1, seq2)])

    return (2 * sim12) / denominator