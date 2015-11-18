from copy import deepcopy
import numpy as np
from analysis import *
from distance import *
from permutation_compstats import *

__all__ = ['vxmatch_siteAnalysis',
           'vxmatch_globalAnalysis',
           'vxmatch_maxt_globalAnalysis',
           'gwj_globalAnalysis']

class vxmatch_siteAnalysis(siteAnalysis):
    """Local sieve analysis method testing if the breakthrough matches the
    vaccine insert at each site (dist is [N x sites])"""

    methodName = 'vxmatch_site'
    def computeDistance(self, params = None):
        """Creates a new distance matrix given the params,
        calling the common _vxmatch function and wiping out previously stored results"""
        self.results.params = deepcopy(params)
        self.results.dist = _vxmatch_distance(self.data.insertSeq, self.data.seqDf, params)

class vxmatch_globalAnalysis(globalAnalysis):
    """Sieve analysis class using a whole-protein distance
    The aggregation across sites happens in prepareDist() which expects
    dist to be [N x sites] but creates fileterdDist as [ N ]
    This is so that the site-wise filters can be applied

    Note: though the distance function is the same as the vxmatch_site analysis, there are other
    specific global functions inherited from globalAnalysis though"""
    
    methodName = 'vxmatch_global'
    def computeDistance(self, params=None):
        """Creates a new distance matrix given the params (wiping out previously stored results)
        This should be identical to the site-wise vxmatch distance functionm since
        filters need to be applied on a site-wise distance matrix (before site aggregation)"""
        self.results.params = deepcopy(params)
        self.results.dist = _vxmatch_distance(self.data.insertSeq, self.data.seqDf, self.results.params)

class gwj_globalAnalysis(globalAnalysis):
    """Global sieve analysis using same distance function as vxmatch_global
    except that the test is a permutation test on the sum of the squared tstatistics at each site
    This is what we call 'Global GWJ' in the RV144 comprehensive sieve manuscript"""
    
    methodName = 'gwj_global'
    comparisonStat = staticmethod(compstat_tstat_sumsquared)
    remoteComparisonStat = 'compstat_tstat_sumsquared'

    def computeDistance(self, params=None):
        """
        Creates a new distance matrix given the params (wiping out previously stored results)
        This should be identical to the site-wise vxmatch distance functionm since
        filters need to be applied on a site-wise distance matrix (before site aggregation)
        """
        self.results.params = deepcopy(params)
        self.results.dist = _vxmatch_distance(self.data.insertSeq, self.data.seqDf, self.results.params)
    def prepareDist(self):
        """Do not compute the summary stat since it needs to be done in the compstat perm test"""
        self.results.filteredDist = deepcopy(self.results.dist.as_matrix())
        self.results.filteredDist[~self.results.distFilter] = np.nan

class vxmatch_maxt_globalAnalysis(sieveAnalysis):
    """Sieve analysis class using a whole-protein distance
    The aggregation across sites happens in prepareDist() which expects
    dist to be [N x sites] but creates fileterdDist as [ N ]
    This is so that the site-wise filters can be applied"""
    
    methodName = 'vxmatch_maxt_global'

    """Here I can use the global scan compstat because filteredDist is 2D [ptid x sites] and we want the maxt over the sites"""
    comparisonStat = staticmethod(compstat_tstat_globalscan)
    remoteComparisonStat = 'compstat_tstat_globalscan'
    def prepareDist(self):
        """Turn the dist DataFrame into a masked ndarray based on the filter
        For a "maxt global" distance leave the distance matrix [N x sites]
        so that filteredDist is also [N x sites]. max() will get applied after group comparison"""
        
        """Use a ndarray instead of masked array to improve speed"""
        self.results.filteredDist = deepcopy(self.results.dist.as_matrix())
        self.results.filteredDist[~self.results.distFilter] = np.nan

    def computeDistance(self, params = None):
        """Creates a new distance matrix given the params (wiping out previously stored results)
        This should be identical to the site-wise vxmatch distance functionm since
        filters need to be applied on a site-wise distance matrix (before site aggregation)"""
        self.results.params = deepcopy(params)
        self.results.dist = _vxmatch_distance(self.data.insertSeq, self.data.seqDf, self.results.params)