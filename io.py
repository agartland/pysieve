"""
Save and load data, analyses and meta-analyses to and from files (db?)

Not sure how I should organize these functions. Should support both file and db backends
Since I'm just going to pickle both should be easy to implement later
Better to have these external functions that can load data and then instantiate neccessary classes
Confusing though because meta.py has its own load/save functions for loading and saving sets of data/analyses,
    but these are only for the simulations i think
"""
__all__ = ['loadSieve',
           'saveSieve']

import pickle
import os
import pandas as pd

def saveSieve(dataPath, obj, dataFn = None, analysisFn = None):
    """Save sieve analysis results and/or data to a file that can be loaded later
       Results and data will be kept in separate files for efficiency if needed.
       Returns the data and results filenames if successful"""
    if dataFn is None:
        dataFn = _getFilename(dataPath, obj.data, 'pkl')
    if analysisFn is None:
        analysisFn = _getFilename(dataPath, obj, 'pkl')

    """If its an analysis object"""
    if hasattr(obj, 'methodName'):
        isAnalysisObj = True
    else:
        isAnalysisObj = False

    if isAnalysisObj:
        analysisClassName = str(obj.__class__).split('.')[-1].replace("'","").replace('>','')
        out = {'methodName':obj.methodName,'analysisClassName':analysisClassName,'results':obj.results}
        with open(analysisFn, 'wb') as fh:
            pickle.dump(out, fh)
    
    """Now save the data"""
    out = {'data':obj.data}
    with open(dataFn, 'wb') as fh:
        pickle.dump(out, fh)
    return dataFn, analysisFn
 
def loadSieve(dataPath, fn, data = None):
    """Load sieve data OR analysis results from a file
    To load data, specify only fn of the data file,
    To load results, specify the pre-loaded data object as data:
        analysisClassObj = loadSieve(DATA_PATH + analysisFn, loadSieve(DATA_PATH + dataFn))

    Parameters
    ----------
    fn : str
        Full path to file
    data : sub-class of pysieve.sieveData
        Specify the data object when loading an analysis object,

    Returns
    -------
    out : sub-class of pysieve.sieveData or pysieve.sieveAnalysis"""

    """Method is compatible across pandas versions and with binary files."""
    out = pd.read_pickle(fn)
   
    """If its an analysis object and we have the data object passed"""
    if 'methodName' in out.keys() and not data is None:
        out['data'] = data
        obj = eval('%s(sievedata = data, sieveresults = results)' % (out['analysisClassName']),globals(),out)
    else:
        obj = out['data']
    return obj

def _getFilename(dataPath, obj, ext):
    """Try to make a filename from as much info is available in the object (data or results)
    Returns the filename"""
    
    """Assume that its an analysis object first"""
    if hasattr(obj,'data') and hasattr(obj.data,'N'):
        isDataObj = False
    else:
        isDataObj = True

    if not isDataObj:
        if obj.data.regionInds is None:
            regStr = 'whole'
        else:
            regStr = '%d_%d' % (obj.data.regionInds[0], obj.data.regionInds[-1])
        
        if obj.results.hlaMethod is None:
            fn = dataPath + '%s/pyresults/%s.%s.%s.%s.%s' % (obj.data.studyName,obj.methodName,obj.data.proteinName,obj.data.insertName,regStr,ext)
        else:
            fn = dataPath + '%s/pyresults/%s.%s.%s.%s.%s.%s' % (obj.data.studyName,obj.methodName,obj.data.proteinName,obj.data.insertName,regStr,obj.results.hlaMethod,ext)
    else:
        """Then its a data object"""
        if obj.regionInds is None:
            regStr = 'whole'
        else:
            regStr = '%d_%d' % (obj.regionInds[0], obj.regionInds[-1])
        if obj.HLAsubset:
            fn = dataPath + '%s/pyresults/data.HLAsubset.%s.%s.%s.%s' % (obj.studyName,obj.proteinName,obj.insertName,regStr,ext)
        else:
            fn = dataPath + '%s/pyresults/data.%s.%s.%s.%s' % (obj.studyName,obj.proteinName,obj.insertName,regStr,ext)
    
    folder,f = os.path.split(fn)

    if not os.path.exists(folder):
        os.mkdir(folder)
    return fn