'''
data.py
contains classes for storing all the data objects for sieve analysis and basic input out routines (e.g. to_csv, to_fasta)

objects include:
    baseData - unsieved set of sequences and HLAs for simulations
    sieveData - minimum dataset needed for sieve analysis: insert, breakthroughs, HLAs, treatment
    simData - a sieveData object containing a simulation property with metadata about the simulation params,date,epitopes etc.
    resultsData - a sieveData object containing results from potentially many sieve analysis methods
    metaResults - contains results of analysis of many sieve datasets
'''

__all__ = ['sieveData',
           'sieveDataMethods']

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped, IUPAC
from Bio.SubsMat.MatrixInfo import blosum90, ident
from StringIO import StringIO

class sieveData(object):
    masterFn = None
    lookupFn = None
    hlaFn = None
    seqFn = None
    mapFn = None

    """lists of unique 2 and 4 digit HLA alleles"""
    uHLA4 = None
    uHLA2 = None
    
    """DataFrame of sequences with index ptid and columns: seq, seqID"""
    seqDf = None
    regionInds = None

    """DataFrame of HLAs with index ptid and columns for all HLA alleles (2 and 4)"""
    hlaDf = None
    hlaFreq = None

    """DataFrame with index ptid and columns: vaccinated, infected, hla"""
    ptidDf = None

    """contains position number as index and hxb2Pos and hxb2aa as columns"""
    mapDf = None

    studyName = None
    proteinName = None
    insertName = None
    
    """sequence strings of the aligned HXB2 and vaccine insert"""
    HXB2 = None
    insertSeq = None

    N = None

    """List of ptids in each group to be used for indexing a df"""
    vacPtid = None
    plaPtid = None
    vacInd = None
    plaInd = None

    temp = {}
    """Signifies to the saving methods that the data may be different than other datasets from the same study"""
    HLAsubset = False

    """Indicates whether the sequence and other site indexed objects have already been sliced by regionInds"""
    isSliced = False

class sieveDataMethods(object):
    data = None
    def __init__(self,sievedata=None):
        if sievedata is None:
            sievedata = sieveData()
        self.data = sievedata

    def isvalidAnalysis(self, proteinName, insertName):
        res = [va for va in s.validAnalyses if va['insertName']==insertName and va['proteinName']==proteinName]
        return len(res) > 0
        
    def to_nexus(self,fn):
        self.to_fasta(fn,fileformat='nexus',sep='_')
    
    def to_fasta(self, fn=None, fileformat='fasta', withHLA=False, withTreatment=False, sep='|', returnString=False):
        """
        >reference|PROTEIN|INSERT
        >ptid|A1|A2|B1|B2 or >ptid|treatment
        >HXB2
        """
        if fn is None:
            fn = '%s.%s.%s.fasta' % (self.data.studyName, self.data.proteinName, self.data.insertName)

        seqRecP = dict(description = '')
        seqP = dict(alphabet = Gapped(IUPAC.protein))

        outList = [SeqRecord(Seq(self.data.insertSeq, **seqP), id = 'reference%s%s%s%s' % (sep,self.data.proteinName,sep,self.data.insertName), **seqRecP),
                   SeqRecord(Seq(self.data.HXB2, **seqP), id = 'HXB2', **seqRecP)]
        tmp = self.data.seqDf.join(self.data.ptidDf)
        for ptid,row in tmp.iterrows():
            treatment = 'vaccine' if row['vaccinated'] else 'placebo'
            idStr = ptid
            if withTreatment:
                idStr += '%s%s' % (sep,treatment)
            if withHLA and 'hla' in row.index and isinstance(row['hla'],basestring):
                idStr += '%s%s' % (sep,sep.join(row['hla']))

            rec = SeqRecord(Seq(row['seq'], **seqP), id = idStr, **seqRecP)
            outList.append(rec)

        if returnString:
            fn = StringIO()
            SeqIO.write(outList, fn, fileformat)
            fn.seek(0)
            return fn.read()
        else:
            SeqIO.write(outList, fn, fileformat)

    def to_treatment_csv(self, fn=None, sep='|', returnString=False):
        if fn is None:
            fn = '%s.%s.%s.trt.csv' % (self.data.studyName, self.data.proteinName, self.data.insertName)

        tmpDf = self.data.seqDf.join(self.data.ptidDf[['vaccinated']], how='left')
        tmpDf['treatment'] = tmpDf.vaccinated.map(lambda s: 'vaccine' if s else 'placebo')
        tmpDf = tmpDf.reset_index()
        tmpDf = tmpDf.rename_axis({'index':'ptid'}, axis=1)
        tmpDf = tmpDf[['ptid','treatment']]
        
        """refPtid = 'reference%s%s%s%s' % (sep,self.data.proteinName,sep,self.data.insertName)
        tmpDf = tmpDf.append({'ptid':refPtid, 'treatment':'reference'}, ignore_index = True)"""
        if returnString:
            fn = StringIO()
            tmpDf.to_csv(fn, index=False)
            fn.seek(0)
            return fn.read()
        else:
            tmpDf.to_csv(fn, index=False)

    def to_mers(self, mersFn=None, nmers=[9], returnList=False):
        allMers = []
        for seq in self.data.seqDf.seq:
            allMers += getMers(seq.replace('-',''), nmers = nmers)
        allMers += getMers(self.data.insertSeq.replace('-',''), nmers = nmers)
        uMers = sorted(list(set(allMers)))
        if returnList:
            return filter(isvalidmer, uMers)
        else:
            with open(mersFn, 'w') as fh:
                for m in uMers:
                    if isvalidmer(m):
                        fh.write('%s\n' % m)
    def to_hla(self, hlaFn = None, returnList = False):
        convert = lambda h: h.replace('_','*')
        if returnList:
            return map(convert,filter(isvalidHLA,self.data.uHLA4))
        else:
            with open(hlaFn,'w') as fh:
                for h in self.data.uHLA4:
                    if isvalidHLA(h):
                        fh.write('%s\n' % convert(h))
    def checkBA(self,ba):
        """Check that all kmers in seqDf and insertSeq are
        present in the binding affinities dict ba, paired with every HLA in hlaDf"""
        tot = 0
        nantot=0

        allMers = []
        for seq in self.data.seqDf.seq:
            allMers += getMers(seq.replace('-',''),nmers=[9])
        allMers += getMers(self.data.insertSeq.replace('-',''),nmers=[9])
        uMers = sorted(list(set(allMers)))
        for m in uMers:
            if isvalidmer(m):
                for h in self.data.uHLA4:
                    if isvalidHLA(h):
                        tot += 1
                        if isnan(ba[(h,m)]):
                            nantot += 1
        print 'Found nan for %d of %d total predictions (%d HLAs, %d mers, %2.0f%% missing)' % (nantot,tot,len(self.data.uHLA4),len(uMers),1e2*nantot/tot)

    def computeDerivedData(self):
        slicestr = lambda yo,ind: ''.join(array([c for c in yo])[array(ind)])

        self.data.N = self.data.seqDf.shape[0]

        """First join ptidDf and seqDf so that plaInd is always a valid boolean index on seqDf"""
        df = self.data.seqDf.join(self.data.ptidDf)
        self.data.vacPtid = df.index[df.vaccinated]
        self.data.plaPtid = df.index[~df.vaccinated]
        """Type of plaInd is ndarray (NOT pd.Series)"""
        self.data.vacInd = df.vaccinated.values.astype(bool)
        self.data.plaInd = (~df.vaccinated).values.astype(bool)

        self.data.ptidDf = df[self.data.ptidDf.columns]
        self.data.seqDf = df[self.data.seqDf.columns]

        """Select region of protein based on regionInds"""
        if not self.data.regionInds is None and not self.data.isSliced:
            rInds = self.data.regionInds
            """Slice seqDf,insertSeq,mapDf,HXB2"""
            for ptid in self.data.seqDf.index:
                seq = self.data.seqDf.seq[ptid]
                self.data.seqDf.seq[ptid] = slicestr(seq,rInds)
            self.data.insertSeq = slicestr(self.data.insertSeq,rInds)

            self.data.mapDf = self.data.mapDf.ix[rInds]
            self.data.mapDf = self.data.mapDf.set_index(arange(len(rInds)))

            self.data.HXB2 = slicestr(self.data.HXB2,rInds)
            self.data.isSliced = True

        """Create df for looking up a site num from HXB2 coordinate"""
        self.data.hxb22site = self.data.mapDf.copy()
        self.data.hxb22site['site'] = self.data.hxb22site.index
        self.data.hxb22site = self.data.hxb22site.set_index('hxb2Pos')
    '''
    TODO: move plotting code to a different file
    def clipXVec(self,hxb2Range = None,vec=None,returnInds=False):
        """Clip seq-axis vector based on an HXB2 coordinate range (eg [70,80])"""
        if hxb2Range is None:
            siteRange = [self.data.mapDf.index[0],self.data.mapDf.index[-1]+1]
        else:
            hxb2Range = [str(c) for c in hxb2Range]
            siteRange = [self.data.mapDf.index[self.data.mapDf.hxb2Pos == hxb2Range[0]],self.data.mapDf.index[self.data.mapDf.hxb2Pos==hxb2Range[1]]+1]
        if returnInds:
            return arange(siteRange[0],siteRange[1])
        else:
            return vec[siteRange[0]:siteRange[1]]
    def plotSeqSpace(self,hxb2Range=None,subst=None,method='tsne',interactive=False,force=False,**kwargs):
        """Plot MDS of sequence space using a substitution matrix. If interactive then returns AnnotationPicker obj"""
        if subst is None:
            subst=blosum90
        seqs=[self.clipXVec(hxb2Range,s) for s in self.data.seqDf.seq]
        df=self.data.ptidDf.join(self.data.seqDf,how='right')
        """uInd has length len(seqs) but indexes into uSeqs"""
        uSeqs,uInd=unique(seqs,return_inverse=True)
        
        group=[]
        for uniqi,s in enumerate(uSeqs):
            tmp=df.vaccinated[uInd==uniqi].unique()
            if len(tmp)==2:
                group.append('both')
            else:
                group.append(tmp[0])
        insertSeq=self.clipXVec(hxb2Range,self.data.insertSeq)
        uSeqs=append(uSeqs,insertSeq)
        group.append('insert')

        recalc=True
        """Recalc if seqMethod doesn't exist or if its different than current method"""
        try:
            if method==self.data.temp['seqMethod']:
                dist=self.data.temp['seqDist']
                xy=self.data.temp['seqXY']
                if xy.shape[0]==len(uSeqs):
                    recalc=False
        except:
            pass

        if recalc or force:
            dist=calcDistanceMatrix(uSeqs,distanceFunc=lambda s1,s2: seq_distance(s1,s2,subst=subst))
            xy=embedDistanceMatrix(dist,method=method)
            self.data.temp['seqDist']=dist
            self.data.temp['seqXY']=xy
            self.data.temp['seqMethod']=method

        freq=objhist(seqs,keys=uSeqs)
        """Make sure the insert has a count of at least 1"""
        if freq[insertSeq]==0:
            freq[insertSeq]=1
        
        if all([f==1 for f in freq.values()]):
            freqVec=[30]*len(freq)
            labels=uSeqs
        else:
            freqVec=scatternorm(array([freq[s] for s in uSeqs]),30,200)
            labels=['%s: %d' % (s,freq[s]) for s in uSeqs]

        if interactive:
            picker=3
        else:
            picker=None
        
        clf()
        scatter(xy[:,0],xy[:,1],s=freqVec,c=[{'insert':'gold','both':'gray',True:'blue',False:'red'}[g] for g in group],picker=picker,**kwargs)
        xticks(())
        yticks(())
        if hxb2Range is None:
            hxb2Range=(self.data.hxb22site.index[0],self.data.hxb22site.index[-1])
        title('MDS Embedding of Sequence space for %s (HXB2 %s-%s)' % (insertSeq,hxb2Range[0],hxb2Range[1]))
        if interactive:
            mp=AnnotationPicker(xy[:, 0], xy[:, 1], labels,weight='bold',color='black',size='x-small')
            return mp

    def plotHLASpace(self,hxb2Range=None,hlaList=None,ba=None,method='tsne',interactive=False,**kwargs):
        """
        Plot an MDS embedding of HLA space
        Original features were nHLAs x nMers
        """
        seqs=[self.clipXVec(hxb2Range,s) for s in self.data.seqDf.seq]
        df=self.data.ptidDf.join(self.data.seqDf,how='right')
        """uInd has length len(seqs) but indexes into uSeqs"""
        uSeqs,uInd=unique(seqs,return_inverse=True)
        
        group=[]
        for uniqi,s in enumerate(uSeqs):
            tmp=df.vaccinated[uInd==uniqi].unique()
            if len(tmp)==2:
                group.append('both')
            else:
                group.append(tmp[0])
        insertSeq=self.clipXVec(hxb2Range,self.data.insertSeq)
        uSeqs=append(uSeqs,insertSeq)
        group.append('insert')

        mers=getMers(insertSeq,nmers=[9])
        dist=empty((len(uSeqs),len(mers)*len(hlaList)))
        for si,s in enumerate(uSeqs):
            for meri, mer in enumerate(getMers(s,nmers=[9])):
                for hlai,h in enumerate(hlaList):
                    pred=ba[(h,mer)]
                    if isnan(pred):
                        pred=15
                    dist[si,int(meri*len(hlaList)+hlai)]=pred    

        xy=embedDistanceMatrix(dist,method=method)
        freq=objhist(seqs,keys=uSeqs)
        """Make sure the insert has a count of at least 1"""
        if freq[insertSeq]==0:
            freq[insertSeq]=1
        
        if all([f==1 for f in freq.values()]):
            freqVec=[30]*len(freq)
            labels=uSeqs
        else:
            freqVec=scatternorm(array([freq[s] for s in uSeqs]),30,200)
            labels=['%s: %d' % (s,freq[s]) for s in uSeqs]

        if interactive:
            picker=3
        else:
            picker=None
        
        clf()
        scatter(xy[:,0],xy[:,1],s=freqVec,c=[{'insert':'gold','both':'gray',True:'blue',False:'red'}[g] for g in group],picker=picker,**kwargs)
        xticks(())
        yticks(())
        if hxb2Range is None:
            hxb2Range=(self.data.hxb22site.index[0],self.data.hxb22site.index[-1])
        title('MDS Embedding of HLA binding space for %s (HXB2 %s-%s)' % (insertSeq,hxb2Range[0],hxb2Range[1]))
        if interactive:
            mp=AnnotationPicker(xy[:, 0], xy[:, 1], labels,weight='bold',color='black',size='x-small')
            return mp
    def plotConservation(self,region=None):
        """Plot entropy/conservation site-wise for vaccine and placebo breakthrough sequences"""
        plotSeqEntropy(self.data.seqDf.seq,region=region)
    '''