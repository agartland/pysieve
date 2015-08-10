from __future__ import division
import numpy as np
from mytstats import tstatistic, nantstatistic

"""
TODO:
    (1) Write then utilize numbized functions in mystats for speedup of permutation testing
"""

__all__=['permfunc',
         'compstat_tstat',
         'compstat_tstat_scan',
         'compstat_tstat_scan_sum',
         'compstat_tstat_globalscan',
         'nonan_compstat_tstat_scan',
         'nonan_compstat_tstat_scan_sum',
         'nonan_compstat_tstat_globalscan',
         'compstat_tstat_sumsquared']

def permfunc(listOfIndices, masked_dist, compstat):
    """Helper function for running a subset of the permutation test (possibly remotely)."""
    """A global dist will be shape [N] while a site dist is [N x sites]
    (EXCEPT for a global_scan which will be [N x nParams] but result is just len=1)"""

    """Run one and look at the dims of the output to pre-allocate results array"""
    indices = listOfIndices[0]
    tmp_result = compstat(masked_dist, **indices)

    if len(masked_dist.shape) > 1 and (not np.isscalar(tmp_result)) and len(tmp_result) > 1:
        results = np.empty((len(listOfIndices), masked_dist.shape[1]))
    else:
        results = np.empty((len(listOfIndices),1))
    for i,indices in enumerate(listOfIndices):
        results[i,:] = compstat(masked_dist, **indices)
    return results

"""
NOTE: aInd and bInd are BOOLEAN indices used to slice the distance matrix
      so the number of ptids in group a is aInd.sum() not aInd.shape[0]
"""

def compstat_tstat(dist, aInd, bInd):
    """
    For local (and global) sieve analysis, compare vaccine and placebo group for each site using a t-statistic
    filteredDist: [ptid x sites] ndarray  (also works for global distances [ptid])

    aInd, bInd: Boolean row index for the two groups
    Returns tstats 1-d ndarray [sites] (for global distances returns a 1d array with the single tstat)
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()
    
    tstat = nantstatistic(a, b, axis = 0, equal_var = False)
    
    """se = np.sqrt((aN-1)*nf.nanvar(a,axis=0)/((aN+bN) - 2) + (bN-1)*nf.nanvar(b,axis=0)/((aN+bN) - 2))
    tstat = (nf.nanmean(a,axis=0) - nf.nanmean(b,axis=0)) / se"""
    if not type(tstat) is np.ndarray:
        tstat = np.array(tstat, ndmin = 1)
    return tstat

def compstat_tstat_scan(dist, aInd, bInd, returnMaxInds = False):
    """
    For local sieve analysis, compare A and B group for each site using a max t-statistic over a parameter space
    filteredDist: [ptid x sites x params] ndarray
    Returns tstat array [sites]

    aInd, bInd: Boolean row index for the two groups
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()

    tstat = nantstatistic(a, b, axis = 0, equal_var = False)
    """se = np.sqrt((aN-1)*nf.nanvar(a,axis=0)/((aN+bN) - 2) + (bN-1)*nf.nanvar(b,axis=0)/((aN+bN) - 2))
    tstat = (nf.nanmean(a,axis=0) - nf.nanmean(b,axis=0)) / se"""

    """tstat.shape --> [sites x params]"""
    sitesNani = np.all(np.isnan(tstat), axis=1)
    """For sites with all nans across params, set all to 0. this makes maxi = 0"""
    tstat[sitesNani,:] = 0
    """Zeros are better than returning nan because if this perm produces a nan
       result then it is not as extreme as observed (which is probably also nan)"""
    maxi = np.nanargmax(np.abs(tstat), axis=1)
    inds = np.ravel_multi_index((np.arange(maxi.shape[0]),maxi), tstat.shape)
    if not returnMaxInds:
        return tstat.flat[inds]
    else:
        return tstat.flat[inds], maxi

def compstat_tstat_scan_sum(dist, aInd, bInd, returnMaxInds = False):
    """Calls compstat_tstat_scan() and then sums across sites.
    Returns the global-fied distance and the optimal params at each site/kmer"""
    tstat,maxi = compstat_tstat_scan(dist, aInd, bInd, returnMaxInds = True)
    """tstats.shape --> [sites]
    maxi.shape --> [sites] with values [0, nParams]"""

    """Shouldn't need nansum here because sites that were nan at all params are defined to be 0
    and all other sites will have at least one non-nan param value"""
    if not returnMaxInds:
        return np.sum(tstat)
    else:
        return np.sum(tstat), maxi

def compstat_tstat_globalscan(dist, aInd, bInd, returnMaxInds = False):
    """
    For global sieve analysis, compare A and B group for each site using a max t-statistic over a parameter space
    filteredDist: [ptid x params] ndarray
    Returns tstat array 1d array with a single element, the tstat

    aInd, bInd: Boolean row index for the two groups

    Also can be used for maxt_global methods where filteredDist is [ptid x sites] and we want to calc a tstat at each site
    and then take the max tstat across sites
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()

    tstat = nantstatistic(a, b, axis = 0, equal_var = False)
    """se = np.sqrt((aN-1)*nf.nanvar(a,axis=0)/((aN+bN) - 2) + (bN-1)*nf.nanvar(b,axis=0)/((aN+bN) - 2))
    tstat = (nf.nanmean(a,axis=0) - nf.nanmean(b,axis=0)) / se"""

    allNan = np.all(np.isnan(tstat), axis = 0)
    """If all nans across params, set to 0. maxi will be 0"""
    if allNan:
        if not returnMaxInds:
            return np.array(np.nan, ndmin=1)
        else:
            return np.array(np.nan, ndmin=1),0
    else:
        """tstat.shape --> [params]"""
        maxi = np.nanargmax(np.abs(tstat), axis=0)
        if not returnMaxInds:
            return np.array(tstat[maxi], ndmin=1)
        else:
            return np.array(tstat[maxi], ndmin=1),maxi

def nonan_compstat_tstat_scan(dist, aInd, bInd, returnMaxInds = False):
    """
    For local sieve analysis, compare A and B group for each site using a max t-statistic over a parameter space
    filteredDist: [ptid x sites x params] ndarray
    Returns tstat array [sites]

    aInd, bInd: Boolean row index for the two groups
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()

    tstat = tstatistic(a, b, axis = 0, equal_var = False)
    """se = np.sqrt((aN-1)*np.var(a,axis=0)/((aN+bN) - 2) + (bN-1)*np.var(b,axis=0)/((aN+bN) - 2))
    tstat = (np.mean(a,axis=0) - np.mean(b,axis=0)) / se"""

    """Even in the nonan cases, the tstat can be nan if there is no variation in either group (divide by zero)"""
    sitesNani = np.all(np.isnan(tstat), axis=1)
    """For sites with all nans across params, set all to 0. this makes maxi = 0"""
    tstat[sitesNani,:] = 0
    """Zeros are better than returning nan because if this perm produces a nan
       result then it is not as extreme as observed (which is probably also nan)"""

    maxi = np.nanargmax(np.abs(tstat), axis=1)
    inds = np.ravel_multi_index((np.arange(maxi.shape[0]), maxi), tstat.shape)
    if not returnMaxInds:
        return tstat.flat[inds]
    else:
        return tstat.flat[inds], maxi

def nonan_compstat_tstat_scan_sum(dist, aInd, bInd, returnMaxInds = False):
    """Calls compstat_tstat_scan() and then sums across sites.
    Returns the global-fied distance and the optimal params at each site/kmer"""
    tstat,maxi = nonan_compstat_tstat_scan(dist, aInd, bInd, returnMaxInds = True)
    """tstats.shape --> [sites]
    maxi.shape --> [sites] with values [0, nParams]"""

    if not returnMaxInds:
        return np.sum(tstat)
    else:
        return np.sum(tstat), maxi

def nonan_compstat_tstat_globalscan(dist, aInd, bInd, returnMaxInds = False):
    """
    For global sieve analysis, compare A and B group for each site using a max t-statistic over a parameter space
    filteredDist: [ptid x params] ndarray
    Returns tstat array 1d array with a single element, the tstat

    aInd, bInd: Boolean row index for the two groups

    Also can be used for maxt_global methods where filteredDist is [ptid x sites] and we want to calc a tstat at each site
    and then take the max tstat across sites
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()

    tstat = tstatistic(a, b, axis = 0, equal_var = False)
    """se = np.sqrt((aN-1)*np.var(a,axis=0)/((aN+bN) - 2) + (bN-1)*np.var(b,axis=0)/((aN+bN) - 2))
    tstat = (np.mean(a,axis=0) - np.mean(b,axis=0)) / se"""

    allNan = np.all(np.isnan(tstat), axis=0)
    """If all nans across params, set to 0. maxi will be 0"""
    if allNan:
        if not returnMaxInds:
            return np.array(np.nan, ndmin=1)
        else:
            return np.array(np.nan, ndmin=1), 0
    else:
        """tstat.shape --> [params]"""
        maxi = np.nanargmax(np.abs(tstat), axis=0)
        if not returnMaxInds:
            return np.array(tstat[maxi], ndmin=1)
        else:
            return np.array(tstat[maxi], ndmin=1), maxi
    
def compstat_tstat_sumsquared(dist, aInd, bInd):
    """
    For global GWJ sieve analysis, compare vaccine and placebo group
    using the sum of the squared t-statistics at each site
    
    filteredDist: [ptid x sites] ndarray

    aInd, bInd: Boolean row index for the two groups
    Returns a 1d array with the single summ (tstat)**2)
    """
    a = dist[aInd]
    b = dist[bInd]

    aN = aInd.sum()
    bN = bInd.sum()

    """Compute site-wise tstatistic"""
    tstat = nantstatistic(a, b, axis = 0, equal_var = False)

    """se = np.sqrt((aN-1)*nf.nanvar(a,axis=0)/((aN+bN) - 2) + (bN-1)*nf.nanvar(b,axis=0)/((aN+bN) - 2))
    tstat = (nf.nanmean(a,axis=0) - nf.nanmean(b,axis=0)) / se"""
    
    sst = nf.nansum(tstat*tstat)
    """If we don't square the tstats then we are looking only for canonical sieve effects. Note also that the test
    becomes one-sided by doing this"""
    #sst = nf.nansum(tstat)

    if not type(sst) is np.ndarray:
        sst = np.array(sst, ndmin=1)
    return sst
