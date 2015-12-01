from __future__ import division
from distance import _nepitope_distance
import numpy as np
from epitope_hotspot import computeHLAHotspots
from analysis_substbased import vxmatch_siteAnalysis

"""
filters.py

Define several filters for sieveData objects to limit the number of sites that are tested
Each filter takes a sieveData object and parameters (also possibly prepared BA matrices)
Returns a boolean array (filter) that is the same dimensions as dist [PTIDs x nSites] or [PTIDs x nSites x nParams]
For a 'filter', False indicates invalid sites/params, while for a mask (used by masked_array) False indicates valid sites/params
    These functions return filters, which allows them to be ANDed together for the intersection of sites/params
"""

__all__ = ['diversityFilter',
           'hotspotFilter',
           'hotspotKmerFilter',
           'binderEscapeKmerFilter',
           'distFilter']

def diversityFilter(sd, minMM = 5, minM = 5, params = {}):
    """Diversity conservation filter based on the insertSeq (False indicates invalid kmers/params)
    Takes a sieveData obj (not the analysis object)"""
    sa = vxmatch_siteAnalysis(sd)
    sa.computeDistance(params = params)

    mm = sa.results.dist
    mmtot = mm.sum(axis = 0)
    mtot = (1 - mm).sum(axis = 0)
    
    mask = ~((mmtot >= minMM) & (mtot >= minM))
    
    reps = np.ones(len(sa.results.dist.shape), dtype = int)
    reps[0] = sa.results.dist.shape[0]

    """Add back the PTID dimension prior to tiling (which shouldn't be neccessary but maybe there's a bug?)"""
    remask = np.array(~mask, ndmin = len(sa.results.dist.shape))
    return np.tile(remask, reps)

def hotspotFilter(sd, hlas, ba, seqs = None, alpha = 0.25, bindingThreshold = 6, nPerms = 1e2):
    """Filter out sites that are not in top fraction of HLA epitope hotspots (False indicates invalid sites/params)"""
    if seqs is None:
        seqs = [sd.insertSeq]
    kmerPdf,kmerp,sitePdf,sitep = computeHLAHotspots(seqs, hlas, ba, bindingThreshold = bindingThreshold, nPerms = nPerms)
    filt = sitep <= alpha
    reps = np.ones(len(sa.results.dist.shape), dtype = int)
    reps[0] = sa.results.dist.shape[0]

    """Add back the PTID dimension prior to tiling (which shouldn't be neccessary but maybe there's a bug?)"""
    remask = np.array(filt, ndmin = len(sa.results.dist.shape))
    return np.tile(remask,reps)

def hotspotKmerFilter(sd, sa, hotspotPct = 0.1):
    """Filter out kmers that are not in top fraction of HLA epitope hotspots (False indicates invalid kmers/params)"""
    binders = _nepitope_distance(sa.results.insertBA,sa.results.btBA,sa.results.params).sum(axis=0)
    ranki = argrank(binders)
    mask = ((ranki/len(binders)) > hotspotPct)
    
    reps = ones(len(sa.results.dist.shape),dtype=int)
    reps[0] = sa.results.dist.shape[0]

    """Add back the PTID dimension prior to tiling (which shouldn't be neccessary but maybe there's a bug?)"""
    remask = array(~mask,ndmin=len(sa.results.dist.shape))
    return tile(remask,reps)

def binderEscapeKmerFilter(sd,sa,binderCriteria=30,escapeCriteria=5):
    """Filter out kmers that don't have the minimum number of predicted insert epitopes and breakthrough escapes (blinded)
    Binder/escape threshold exclude kmers < threshold
    Returns a filter with same dimensions as dist (False indicates invalid kmers/params)"""

    """binders is [PTIDs x nSites] or for  a param scan distance [PTIDs x nSites x nParams]"""
    binders = _nepitope_distance(sa.results.insertBA, sa.results.btBA, sa.results.params).sum(axis=0)
    escapes = sa.results.dist.sum(axis=0)
    mask = (binders < binderCriteria) | (escapes < escapeCriteria)
    
    """Add back the PTID dimension prior to tiling (which shouldn't be neccessary but maybe there's a bug?)"""
    remask = np.array(~mask, ndmin = len(sa.results.dist.shape))
    
    """Create a reps arg for tiling with nPTIDs in first dim and ones in the remaining dims"""
    reps = np.ones(len(sa.results.dist.shape), dtype = int)
    reps[0] = sa.results.dist.shape[0]
    out = np.tile(remask, reps)
    return out

def distFilter(sd, sa, distCriteria = 5):
    """Filter on blinded distane alone.
    Returns a filter with same dimensions as dist (False indicates invalid kmers/params)"""
    
    dist = sa.results.dist.sum(axis=0)
    mask = dist < distCriteria
    
    """Add back the PTID dimension prior to tiling (which shouldn't be neccessary but maybe there's a bug?)"""
    remask = np.array(~mask, ndmin = len(sa.results.dist.shape))
    
    """Create a reps arg for tiling with nPTIDs in first dim and ones in the remaining dims"""
    reps = np.ones(len(sa.results.dist.shape), dtype = int)
    reps[0] = sa.results.dist.shape[0]
    out = np.tile(remask, reps)
    return out
