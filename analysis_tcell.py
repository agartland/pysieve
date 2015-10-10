from copy import deepcopy
import numpy as np
import itertools
import pandas as pd
from decimal import Decimal

from analysis import *
from distance import *
from permutation_compstats import *

__all__ = ['relative_binding_escape_kmerAnalysis',
           'binding_escape_kmerAnalysis',
           'indel_escape_kmerAnalysis',
           'indel_binding_escape_kmerAnalysis',
           'binding_scan_kmerAnalysis',
           'indel_binding_scan_kmerAnalysis',
           'relative_binding_escape_globalAnalysis',
           'binding_escape_globalAnalysis',
           'indel_escape_globalAnalysis',
           'indel_binding_escape_globalAnalysis',
           'binding_scan_globalAnalysis',
           'indel_binding_scan_globalAnalysis',
           'binding_kmerscan_globalAnalysis']

class relative_binding_escape_kmerAnalysis(siteAnalysis):
    """Similar to binding escape except that it uses the median of non-binding HLAs as an escape threshold
    Performed on the kmer level (the sites dimension is kmer start positions)"""
    methodName = 'relative_binding_escape_kmer'

    def initialize(self, ba, params = {'nmer':9}):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.params['fullBT'] = True
        self.results.analysisMethod = self.methodName
        
        out = _prepareBA(self.data, ba, self.results.params)
        self.results.insertBA, self.results.btBA, self.results.hlaMethod = out
    def computeDistance(self, params={'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.params['fullBT'] = True
        self.results.dist = _relative_binding_escape_distance(self.data.insertSeq,
                                    self.results.insertBA,
                                    self.data.seqDf,
                                    self.results.btBA,
                                    self.data.hlaDf[self.data.uHLA4].ix[self.data.seqDf.index],
                                    self.results.params, None)

class binding_escape_kmerAnalysis(siteAnalysis):
    """'HLA binding escape count' performed on the kmer level
    (the sites dimension is kmer start positions)"""

    methodName='binding_escape_kmer'

    def initialize(self, ba, params = {'nmer':9}):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.analysisMethod = self.methodName
        out = _prepareBA(self.data, ba, self.results.params)
        self.results.insertBA, self.results.btBA,self.results.hlaMethod = out
    def computeDistance(self, params = {'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _binding_escape_distance(self.data.insertSeq,
                                                     self.results.insertBA,
                                                     self.data.seqDf,
                                                     self.results.btBA,
                                                     params)

class indel_escape_kmerAnalysis(binding_escape_kmerAnalysis):
    """HLA indel escape count' performed on the kmer level"""
    methodName = 'indel_escape_kmer'

    def computeDistance(self, params = {'binding':5,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _indel_escape_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   params)

class indel_binding_escape_kmerAnalysis(binding_escape_kmerAnalysis):
    """HLA indel + binding escape count' performed on the kmer level"""
    methodName = 'indel_binding_escape_kmer'

    def computeDistance(self, params = {'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _indel_escape_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   params)
        self.results.dist += _binding_escape_distance(self.data.insertSeq,
                                                      self.results.insertBA,
                                                      self.data.seqDf,
                                                      self.results.btBA,
                                                      params)

class binding_scan_kmerAnalysis(siteScanAnalysis):
    """
    Similar to the HLA binding escape count distance
    except that a parameter scan is performed (independently for each kmer)
    over the binding and escape threshold parameter ranges
    Default ranges are specified here in computeDistance
    """
    methodName = 'binding_scan_kmer'

    @property
    def paramPairs(self):
        bindingRange = np.arange(self.results.params['bindingStart'], self.results.params['bindingEnd'], self.results.params['gridDt'])
        escapeRange = np.arange(self.results.params['escapeStart'], self.results.params['escapeEnd'], self.results.params['gridDt'])
        return [(Decimal('%1.1f' % b), Decimal('%1.1f' % e)) for b,e in itertools.product(bindingRange, escapeRange)]

    def initialize(self, ba, params):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.analysisMethod = self.methodName
        out = _prepareBA(self.data, ba, params)
        self.results.insertBA, self.results.btBA, self.results.hlaMethod = out

    def computeDistance(self, params):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _binding_scan_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   self.paramPairs,
                                                   params['minDelta'],
                                                   params['nmer'])

    def computeObserved(self, distFilter=None):
        """
        Computes the observed distance between treatment groups using the filteredDist and the comparisonStat
        For a parameter scan analysis also creates in params a scannedBinding and scannedEscape which are the parmater pairs
        that were optimal for each site
        """
        if distFilter is None:
            """If no filter exists then create an all inclusive filter"""
            distFilter = np.ones(self.results.dist.shape, dtype = bool)
        """Filter out nan distances (per PTID, not the whole site)"""
        """When dist is an ndarray, don't try to use the .values attribute"""
        distFilter[np.isnan(self.results.dist)] = False
        self.results.distFilter = distFilter
        self.prepareDist()
        self.results.observed,maxi = self.comparisonStat(self.results.filteredDist,
                                                         aInd = self.data.vacInd,
                                                         bInd = self.data.plaInd,
                                                         returnMaxInds = True)
        paramPairs = self.paramPairs
        distShape = self.results.dist.shape

        scannedDist = np.nan * np.zeros(distShape[:2])
        for i,mxi in enumerate(maxi):
            scannedDist[:,i] = self.results.dist[:,i,mxi]
        self.results.scannedDist = pd.DataFrame(scannedDist,
                                                index = self.data.seqDf.index,
                                                columns = np.arange(distShape[1]))
        self.results.params['scannedBinding'] = np.array([paramPairs[i][0] for i in maxi])
        self.results.params['scannedEscape'] = np.array([paramPairs[i][1] for i in maxi])

class indel_binding_scan_kmerAnalysis(binding_scan_kmerAnalysis):
    """HLA indel + binding escape count performed on the kmer level with a parameter scan"""
    methodName = 'indel_binding_scan_kmer'

    def computeDistance(self,params):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _indel_scan_distance(self.data.insertSeq,
                                                 self.results.insertBA,
                                                 self.data.seqDf,
                                                 self.results.btBA,
                                                 params['nmer'],
                                                 self.paramPairs,
                                                 params['minDelta'])
        self.results.dist += _binding_scan_distance(self.data.insertSeq,
                                                    self.results.insertBA,
                                                    self.data.seqDf,
                                                    self.results.btBA,
                                                    self.paramPairs,
                                                    params['minDelta'],
                                                    params['nmer'])

class relative_binding_escape_globalAnalysis(globalAnalysis):
    """Global sieve analysis using the relative binding escape count distance (i.e. Allan's method)"""
    methodName = 'relative_binding_escape_global'

    def initialize(self, ba, params = {'nmer':9,'fullBT':True}):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.params['fullBT'] = True
        self.results.analysisMethod = self.methodName
        out = _prepareBA(self.data, ba, self.results.params)
        self.results.insertBA, self.results.btBA, self.results.hlaMethod = out

    def computeDistance(self, params = {'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        if self.results.params is None:
            self.results.params = deepcopy(params)
        self.results.params['fullBT'] = True
        self.results.dist = _relative_binding_escape_distance(self.data.insertSeq,
                                        self.results.insertBA,
                                        self.data.seqDf,
                                        self.results.btBA,
                                        self.data.hlaDf[self.data.uHLA4].ix[self.data.seqDf.index],
                                        self.results.params,None)
class binding_escape_globalAnalysis(globalAnalysis):
    """Global sieve analysis using the escape count distance"""
    methodName = 'binding_escape_global'

    def initialize(self,ba,params={'nmer':9}):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.analysisMethod = self.methodName
        out = _prepareBA(self.data, ba, self.results.params)
        self.results.insertBA, self.results.btBA, self.results.hlaMethod = out

    def computeDistance(self, params = {'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _binding_escape_distance(self.data.insertSeq,
                                                     self.results.insertBA,
                                                     self.data.seqDf,
                                                     self.results.btBA,
                                                     self.results.params)

class indel_escape_globalAnalysis(binding_escape_globalAnalysis):
    """Global sieve analysis using the indel escape count distance"""
    methodName = 'indel_escape_global'

    def computeDistance(self, params = {'binding':5,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _indel_escape_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   self.results.params)

class indel_binding_escape_globalAnalysis(binding_escape_globalAnalysis):
    """Global sieve analysis using the indel + binding escape count distance"""
    methodName = 'indel_binding_escape_global'

    def computeDistance(self, params = {'binding':5,'escape':7,'nmer':9}):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _indel_escape_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   self.results.params)
        self.results.dist += _binding_escape_distance(self.data.insertSeq,
                                                      self.results.insertBA,
                                                      self.data.seqDf,
                                                      self.results.btBA,
                                                      self.results.params)

class binding_scan_globalAnalysis(globalScanAnalysis):
    """Global sieve analysis using the escape count distance with a parameter scan"""
    methodName = 'binding_scan_global'

    def initialize(self, ba, params = {'nmer':9}):
        """Calls the private-like function defined in analysis.py"""
        self.dropMissingHLA()
        self.results.params = deepcopy(params)
        self.results.analysisMethod = self.methodName
        out = _prepareBA(self.data, ba, self.results.params)
        self.results.insertBA, self.results.btBA, self.results.hlaMethod = out

    @property
    def paramPairs(self):
        bindingRange = np.arange(self.results.params['bindingStart'], self.results.params['bindingEnd'], self.results.params['gridDt'])
        escapeRange = np.arange(self.results.params['escapeStart'], self.results.params['escapeEnd'], self.results.params['gridDt'])
        return [(Decimal('%1.1f' % b),Decimal('%1.1f' % e)) for b,e in itertools.product(bindingRange, escapeRange)]

    def computeDistance(self,params):
        """Creates a new distance matrix given the params (wiping out previously stored results)"""
        self.results.params = deepcopy(params)
        self.results.dist = _binding_scan_distance(self.data.insertSeq,
                                                   self.results.insertBA,
                                                   self.data.seqDf,
                                                   self.results.btBA,
                                                   self.paramPairs,
                                                   params['minDelta'],
                                                   params['nmer'])

    def computeObserved(self, distFilter = None):
        """Computes the observed distance between treatment groups using the filteredDist and the comparisonStat
        For a parameter scan analysis also creates in params a scannedBinding and scannedEscape which are the parmater pairs
        that were optimal for each site"""
        if distFilter is None:
            """If no filter exists then create an all inclusive filter"""
            distFilter = np.ones(self.results.dist.shape, dtype = bool)
        """Filter out nan distances (per PTID, not the whole site)"""
        """When dist is an ndarray, don't try to use the .values attribute"""
        distFilter[np.isnan(self.results.dist)] = False
        self.results.distFilter = distFilter
        self.prepareDist()
        self.results.observed, maxi = self.comparisonStat(self.results.filteredDist,
                                                          aInd = self.data.vacInd,
                                                          bInd = self.data.plaInd,
                                                          returnMaxInds = True)
        optimalParams = self.paramPairs[maxi]
        distShape = self.results.dist.shape

        """dist and scannedDist are [N x sites] but scannedDist is missing params, which is good for plotting optimal later"""
        scannedDist = self.results.dist[:,:,maxi]
        self.results.scannedDist = pd.DataFrame(scannedDist,
                                                index = self.data.seqDf.index,
                                                columns = np.arange(distShape[1]))
        self.results.params['scannedBinding'] = optimalParams[0]
        self.results.params['scannedEscape'] = optimalParams[1]

class indel_binding_scan_globalAnalysis(binding_scan_globalAnalysis):
    """Global sieve analysis using the indel + binding escape count distance"""
    methodName = 'indel_binding_scan_global'

    def computeDistance(self, params):
        self.results.params = deepcopy(params)
        self.results.dist = _indel_scan_distance(self.data.insertSeq,
                                                 self.results.insertBA,
                                                 self.data.seqDf,
                                                 self.results.btBA,
                                                 params['nmer'],
                                                 self.paramPairs,
                                                 params['minDelta'])
        self.results.dist += _binding_scan_distance(self.data.insertSeq,
                                                    self.results.insertBA,
                                                    self.data.seqDf,
                                                    self.results.btBA,
                                                    self.paramPairs,
                                                    params['minDelta'],
                                                    params['nmer'])

class binding_kmerscan_globalAnalysis(binding_scan_globalAnalysis):
    """Very similar to the local binding_scan_kmer analysis except the distance is summed across kmers
    Differs from the binding_scan_global method because the optimal thresholds are chosen for each kmer
        instead of one pair of thresholds for the whole protein (this make it slower and probably less powerful)

    Similar to the HLA binding escape count distance
    except that a parameter scan is performed (independently for each kmer)
    over the binding and escape threshold parameter ranges
    Default ranges are specified her in computeDistance

    By using a NONAN compstat all nans must be removed in filteredDist by computeObserved()"""
    methodName = 'binding_kmerscan_global'
    comparisonStat = staticmethod(nonan_compstat_tstat_scan_sum)
    remoteComparisonStat = 'nonan_compstat_tstat_scan_sum'

    def prepareDist(self):
        """Turn the dist DataFrame into a nan-masked ndarray based on the filter
        
        For a "global" distance this also mean taking the distance matrix [N x sites x params]
        and computing a summary statistc so that filteredDist is [N x params]

        BUT for a kmerscan global method the site-aggregation will not happen until the parameter optimization
        so the summing doesn't happen here (hence this different prepareDist())

        The only difference relative to the prepareDist() in globalScanAnalysis is the sum() line that has been commented out"""

        """Since the filter will be of shape [PTIDS x nSites x nparams] use it as is"""
        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist)
        self.results.filteredDist[~self.results.distFilter] = np.nan

        """Removing nans is fine for escape count-based distances and should speed up computation when paired with NONAN compstats"""
        self.results.filteredDist[np.isnan(self.results.filteredDist)] = 0

        #self.results.filteredDist=nansum(self.results.filteredDist,axis=1)

    def computeObserved(self, distFilter = None):
        """
        Computes the observed distance between treatment groups using the filteredDist and the comparisonStat
        For a parameter scan analysis also creates in params a scannedBinding and scannedEscape which are the parmater pairs
        that were optimal for each site
        """
        if distFilter is None:
            """If no filter exists then create an all inclusive filter"""
            distFilter = np.ones(self.results.dist.shape, dtype = bool)
        """Filter out nan distances (per PTID, not the whole site)"""
        """When dist is an ndarray, don't try to use the .values attribute"""
        distFilter[np.isnan(self.results.dist)] = False
        self.results.distFilter = distFilter
        self.prepareDist()
        
        self.results.observed,maxi = self.comparisonStat(self.results.filteredDist,
                                                         aInd = self.data.vacInd,
                                                         bInd = self.data.plaInd,
                                                         returnMaxInds = True)
        paramPairs = self.paramPairs
        distShape = self.results.dist.shape

        """scannedDist is [ptid x sites] (same shape as that for binding_scan_globalAnalysis)"""
        scannedDist = np.nan * np.zeros(distShape[:2])
        for i,mxi in enumerate(maxi):
            scannedDist[:,i] = self.results.dist[:,i,mxi]
        self.results.scannedDist=pd.DataFrame(scannedDist,
                                              index = self.data.seqDf.index,
                                              columns = np.arange(distShape[1]))
        self.results.params['scannedBinding'] = np.array([paramPairs[i][0] for i in maxi])
        self.results.params['scannedEscape'] = np.array([paramPairs[i][1] for i in maxi])