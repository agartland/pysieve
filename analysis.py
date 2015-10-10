import pandas as pd
import numpy as np
from copy import deepcopy
import logging
import statsmodels.api as sm
import os.path as op

from permutation_compstats import *

__all__=['sieveResults',
         'sieveAnalysis',
         'siteAnalysis',
         'siteScanAnalysis',
         'globalAnalysis',
         'globalScanAnalysis']

class sieveResults(object):
    """Container for sieve results. Each method will have its own results object
    (because it would be rare to need to compare across methods and it will be more efficient
    to have the results split up if they are going to be passed around the cluster)
    """

    """Each of these can be shaped to contain X for each ptid and or each site/kmer"""
    dist = None
    params = None
    observed = None
    permutations = None
    pvalue = None

    """distFilter is False for invalid sites/params (opposite of that required for mask in masked_array)"""
    distFilter = None
    filteredDist = None
    unfilterFunc = None

    nPerms = 1e2
    randomSeed = 22

    hlaMethod = None
    analysisMethod = None

    temp = {}
    

class sieveAnalysis(object):
    def __init__(self, sievedata, sieveresults = None):
        self.data = sievedata
        if sieveresults is None:
            self.results = sieveResults()
        else:
            """Or consider ability to load from a file given the file location"""
            self.results = sieveresults

        """Set up logger for printing progress of long-running simulations and analysis"""
        self.logger = logging.getLogger('pysieve.analysis')
        self.logger.setLevel(logging.INFO)

        formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(message)s')

        """Create console handler and set level to debug"""
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
    def initialize(self, *args, **kwargs):
        """Used to prepare insertBA and btBA or any other computation prior to computing the distance itself
        Should be called after instantiation"""
        self.results.analysisMethod = self.methodName
    def dropMissingHLA(self):
        """Drops ptids from hlaDf and ptidDf and seqDf if they are missing any HLA alleles"""
        keep = self.data.ptidDf.index[~self.data.ptidDf.hla.isnull()]
        nbefore = self.data.ptidDf.shape[0]
        self.data.ptidDf = self.data.ptidDf.loc[keep]
        self.data.hlaDf = self.data.hlaDf.loc[keep]
        self.data.seqDf = self.data.seqDf.loc[keep]
        self.data.N = len(keep)
        self.data.vacInd = self.data.ptidDf.vaccinated.values.astype(bool)
        self.data.plaInd = (~self.data.ptidDf.vaccinated).values.astype(bool)
        self.data.vacPtid = self.data.ptidDf.index[self.data.ptidDf.vaccinated]
        self.data.plaPtid = self.data.ptidDf.index[~self.data.ptidDf.vaccinated]
        self.data.HLAsubset = True
        self.logger.info('Dropped %d participants that were missing HLA data (N = %d ---> %d)', nbefore-len(keep), nbefore, len(keep))
    
    def identifierStr(self):
        """Return identifier string (e.g. "RV144 env (MN): vxmatch_site")"""
        return '%s %s (%s): %s' % (self.data.studyName, self.data.proteinName, self.data.insertName, self.results.analysisMethod)
    def thinPermutations(self, n = 1000):
        if not self.results.permutations is None:
            if self.results.permutations.shape[0] > n:
                self.results.permutations = self.results.permutations[:n]
    def to_distance_csv(self, fn = None):
        """Save the distances to a CSV file"""
        if fn is None:
            fn = '%s.%s.%s.%s.distance.csv' % (self.data.studyName, self.data.proteinName, self.data.insertName, self.results.analysisMethod)
        self.to_distance_df().to_csv(fn, index = False, na_rep='NAN')
    def to_distance_df(self):
        """Return distances in a pd.DataFrame()"""
        columns = ['ptid',
                   'display_position',
                   'start_position',
                   'distance_method',
                   'distance']
        try:
            dist = self.results.scannedDist.values
        except:
            dist = self.results.filteredDist

        d = {col:[] for col in columns}

        for ptidi in range(dist.shape[0]):
            for sitei in range(dist.shape[1]):
                d['ptid'].append(self.data.seqDf.index[ptidi])
                d['display_position'].append(self.data.mapDf.hxb2Pos[sitei])
                d['start_position'].append(sitei)
                d['distance_method'].append(self.results.analysisMethod)
                d['distance'].append(dist[ptidi,sitei])

        df = pd.DataFrame(d)[columns]
        return df
        
    def to_results_csv(self, fn=None):
        self.to_csv(fn)
    def to_csv(self, fn=None, fdir=None):
        """Save the results to a CSV file"""
        if fn is None:
            fn = '%s.%s.%s.%s.results.csv' % (self.data.studyName, self.data.proteinName, self.data.insertName, self.results.analysisMethod)
        if not fdir is None:
            fullFn = op.join(fdir, fn)
        else:
            fullFn = fn
        self.to_df().to_csv(fullFn, index  = False, na_rep='NAN')
    def to_df(self):
        """Return results in a pd.DataFrame()"""
        if isinstance(self, globalAnalysis):
            columns = ['distance_method',
                       'mean_placebo_distance',
                       'mean_vaccine_distance',
                       'observed_statistic',
                       'pvalue']
            df = pd.DataFrame(np.zeros((1,len(columns)), dtype = object),index = None, columns = columns)
            try: 
                """Sum across sites"""
                dist = self.results.scannedDist.sum(axis=1)
            except:
                dist = self.results.filteredDist
            if self.results.analysisMethod == 'vxmatch_maxt_global':
                dist = dist.max(axis=1)
            plaDist = dist[self.data.plaInd].mean(axis=0)
            vacDist = dist[self.data.vacInd].mean(axis=0)
            observed = self.results.observed
            if not np.isscalar(observed):
                observed = np.float(observed)
            pvalue = self.results.pvalue
            if not np.isscalar(pvalue):
                pvalue = np.float(pvalue)
            df.iloc[0] = pd.Series([self.results.analysisMethod, plaDist,vacDist,observed,pvalue],index = columns)
        else:
            columns = ['distance_method',
                       'display_position',
                       'start_position',
                       'placebo_distance',
                       'vaccine_distance',
                       'sieve_statistic',
                       'pvalue',
                       'qvalue']
            try:
                dist = self.results.scannedDist.values
            except:
                dist = self.results.filteredDist
            plaDist = dist[self.data.plaInd,:].mean(axis=0)
            vacDist = dist[self.data.vacInd,:].mean(axis=0)
            observed = self.results.observed
            pvalue = self.results.pvalue
            _, qvalue, _, _ = sm.stats.multipletests(pvalue, method = 'fdr_bh', alpha = 0.2)
            
            df = pd.DataFrame(np.zeros((len(plaDist), len(columns)), dtype = object),index = None, columns = columns)
            for sitei,(pdt,vdt,obs,p,q) in enumerate(zip(plaDist,vacDist,observed,pvalue,qvalue)):
                df.iloc[sitei] = pd.Series([self.results.analysisMethod,
                                            self.data.mapDf.hxb2Pos[sitei],
                                            sitei,
                                            pdt,
                                            vdt,
                                            obs,
                                            p,
                                            q],index = columns)
        return df
    def computeDistance(self, params = None):
        pass
    @staticmethod
    def comparisonStat(masked_dist, aInd, bInd):
        pass
    def computeObserved(self, distFilter = None):
        """Compute the observed value of the difference between the vaccine and placebo groups (possibly, for each site)
           Creates filteredDist which is needed to call permutationTest()"""
        if distFilter is None:
            """If no filter exists then create a filter that excludes nans"""
            distFilter = np.ones(self.results.dist.shape, dtype = bool)
        """Filter out nan distances (per PTID, not the whole site)"""
        distFilter[np.isnan(self.results.dist.values)] = False
        self.results.distFilter = distFilter
        """Prepare dist applies the distFilter and creates filteredDist (which also aggregates across sites for global distances)"""
        self.prepareDist()
        tmpFunc = self.comparisonStat
        self.results.observed = tmpFunc(self.results.filteredDist, aInd = self.data.vacInd, bInd = self.data.plaInd)
    def permutationTest(self, nperms = None, clusterClient = None):
        """Permutation test using the comparison stat over the two treatment groups based on distances in self.filteredDist
        dist: can be [ptid x sites] or [ptid] or whatever as long as the comparisonStat can handle filteredDist correctly
        observed: tstat for each site (pd.Series with arange(nSites as index))

        Note: must be preceded by calls to computeDistance() and computeObserved() to set up self.filteredDist"""
        if nperms is None:
            nperms = self.nPerms
        else:
            self.nPerms = nperms

        vacN = self.data.vacInd.sum()
        plaN = self.data.plaInd.sum()

        sitesN = self.results.dist.shape[1]

        """Predetermine the permutation indices"""
        np.random.seed(self.results.randomSeed)
        randInds = []
        for permi in range(nperms):
            ind = np.random.permutation(self.data.N)
            randInds.append({'aInd':ind[:vacN],'bInd':ind[vacN:]})
        
        if clusterClient is None:
            self.results.permutations = permfunc(randInds,self.results.filteredDist,self.comparisonStat)
        else:
            clusterClient[:]['remote_dist'] = self.results.filteredDist
            clusterClient[:].scatter('remote_inds', randInds, block=True)
            clusterClient[:].execute('import pysieve.permutation_compstats', block=True)
            clusterClient[:].execute('remote_result = pysieve.permutation_compstats.permfunc(remote_inds,remote_dist,pysieve.permutation_compstats.%s)' % self.remoteComparisonStat, block=True)
            self.results.permutations = clusterClient[:].gather('remote_result', block=True)

            #self.results.permutations=ParallelFunction(clusterClient[:],permTestFunc,block=True,chunksize=nperms/(2*len(clusterClient)))(randInds)
            
    def computePvalues(self):
        observed = np.tile(np.array(self.results.observed), (self.results.permutations.shape[0],1))
        self.results.pvalue = 1 - (np.abs(self.results.permutations) < np.abs(observed)).sum(axis=0).astype(np.float64)/observed.shape[0]
        self.results.pvalue[self.results.pvalue < (np.float64(1)/np.float64(observed.shape[0]))] = np.float64(1)/observed.shape[0]
    def computeFST(self, nperms = 10, subst = None):
        """Compute a statistic that compares pairwise diversity (PD) within group to PD between groups,
        and performs a permutation test

        TODO: Test since switched to use seqdistance"""
        if subst is None:
            self.logger.info('Using Nan gapScores (i.e. ignoring gap comparisons)')
            subst = seqdistance.matrices.addGapScores(seqdistance.matrices.binarySubst, seqdistance.matrices.binGapScores)

        PDFunc = lambda a,b: np.nanmean(seqdistance.distance_df(a.tolist(),b.tolist(), args = (subst), symteric = True))

        def fst(align, vInd, pInd):
            VP = PDFunc(align.loc[vInd], align.loc[pInd])
            VV = PDFunc(align.loc[vInd], align.loc[vInd])
            PP = PDFunc(align.loc[pInd], align.loc[pInd])
            return (VP - (VV + PP) / 2) / VP

        vacN = self.data.vacInd.sum()
        plaN = self.data.plaInd.sum()
        
        obs = fst(self.data.seqDf.seq, self.data.vacInd, self.data.plaInd)
        
        samps = np.zeros(nperms)
        np.random.seed(self.results.randomSeed)
        for permi in range(nperms):
            ind = np.random.permutation(self.data.N)
            samps[permi] = fst(self.data.seqDf.seq, ind[:vacN], ind[vacN:])
        pvalue = (samps > obs).sum() / nperms
        return obs, pvalue

class siteAnalysis(sieveAnalysis):
    """Base class for all site- or kmer-based analyses (distance matrix has dims [N x nSites])"""
    methodName = 'site'
    comparisonStat = staticmethod(compstat_tstat)
    remoteComparisonStat = 'compstat_tstat'
    def prepareDist(self):
        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist.as_matrix())
        self.results.filteredDist[~self.results.distFilter] = np.nan
        '''
        """Turn the dist DataFrame into a masked ndarray based on the filter"""
        mask=~self.results.distFilter
        self.results.filteredDist=ma.masked_array(data=self.results.dist.as_matrix(),
                                                  mask=mask,
                                                  fill_value=nan)
        '''
    def rankFeatures(self):
        """Use supervised ML methods to identify and rank site-based features that are
        predictive of treatment assignment
        Implemented this as separate analysis methods in ml.py"""
        pass

class siteScanAnalysis(siteAnalysis):
    """Base class for all site- or kmer-based analyses with a parameter scan (distance matrix has dims [N x nSites x nParams])"""
    methodName = 'site_scan'
    comparisonStat = staticmethod(compstat_tstat_scan)
    remoteComparisonStat = 'compstat_tstat_scan'

    def prepareDist(self):
        """Turn the dist DataFrame into a masked ndarray based on the filter"""
        """Since the filter will be of shape [PTIDS x nSites] use it as is"""
        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist)
        self.results.filteredDist[~self.results.distFilter] = np.nan
        '''
        mask=~self.results.distFilter
        self.results.filteredDist=ma.masked_array(data=self.results.dist,mask=mask,fill_value=nan)
        '''
    
class globalAnalysis(sieveAnalysis):
    """Base class for all global analyses (distance matrix has dims [ N ])"""
    methodName = 'global'
    comparisonStat = staticmethod(compstat_tstat)
    remoteComparisonStat = 'compstat_tstat'
    def prepareDist(self):
        """Turn the dist DataFrame into a masked ndarray based on the filter
        For a "global" distance this also mean taking the distance matrix [N x sites]
        and computing a summary statistc so that filteredDist is [N]
        The inherited/default summary statistic is a sum(axis = 1) across sites, divide by the number of valid sites
            resulting in a distance per site scaling"""

        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist.as_matrix())
        self.results.filteredDist[~self.results.distFilter] = np.nan
        self.results.filteredDist = np.nansum(self.results.filteredDist, axis=1) / (~np.isnan(self.results.filteredDist)).sum(axis=1)
        '''
        mask=~self.results.distFilter
        self.results.filteredDist=ma.masked_array(data=self.results.dist.as_matrix(),
                                                  mask=mask,
                                                  fill_value=nan)
        self.results.filteredDist=self.results.filteredDist.sum(axis=1)
        '''


class globalScanAnalysis(globalAnalysis):
    """Base class for all global analyses with a parameter scan
    The distance matrix has dims [N x nSites x nParams]
    But the filteredDist (on which the tstat performs) has dims [N x nParams]"""
    methodName = 'global_scan'
    comparisonStat = staticmethod(nonan_compstat_tstat_globalscan)
    remoteComparisonStat = 'nonan_compstat_tstat_globalscan'

    def prepareDist(self):
        """Turn the dist DataFrame into a masked ndarray based on the filter
        
        For a "global" distance this also mean taking the distance matrix [N x sites x params]
        and computing a summary statistc so that filteredDist is [N x params]
        The inherited/default summary statistic is a sum(axis = 1) across sites (NOT CURRENTLY NORMALIZED BY THE NUMBER OF SITES)"""

        """Since the filter will be of shape [PTIDS x nSites x nparams] use it as is"""
        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist)
        self.results.filteredDist[~self.results.distFilter] = np.nan
        self.results.filteredDist = np.nansum(self.results.filteredDist, axis=1)

        """Removing nans is fine for escape count-based distances and should speed up computation when paired with NONAN compstats"""
        self.results.filteredDist[np.isnan(self.results.filteredDist)] = 0
        
        '''
        mask=~self.results.distFilter
        self.results.filteredDist=ma.masked_array(data=self.results.dist,mask=mask,fill_value=nan)
        self.results.filteredDist=self.results.filteredDist.sum(axis=1)
        '''

        """
        if (isnan(self.results.dist) & (~(mask))).sum()>0:
            raise Exception('Some Nans not masked!')
        a=self.results.filteredDist.sum(axis=1)
        b=nansum(self.results.filteredDist,axis=1)
        if not all(isnan(a)==isnan(b)):
            raise Exception('Sum and nansum dont have same nans')
        a=a[~isnan(a)]
        b=b[~isnan(b)]
        if not all(a==b):
            raise Exception('Sum and nansum are different in values')
        """






