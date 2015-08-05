from __future__ import division
"""
meta.py
contains classes for meta-analysis of multiple simulation datasets

inputs: multiple sieve results objects

methods:
    (1) plotting of sensitivity, specificity, performance etc. by parameter values
    (2) can be distance specific (e.g. make plots about epitopes identified etc.)
"""

import pickle
import sys
import os
import time
from copy import deepcopy
import pandas as pd
import numpy as np
import itertools

from simulation import *
#from evaluation import *

"""
TODO:
    (1) Use logger module instead of print statements
    (2) Work on plotting methods (currently untested and non-functioning
"""

__all__=['simulationMeta',
         'analysisMeta',
         'siteMeta',
         'loadGlobalPvalues',
         'plotPRByP',
         'plotFPR',
         'plotPRatP']

class analysisMeta(object):
    """Class for organizing, saving/loading analyses (single method) of a simulationMeta object"""
    def __init__(self, dataPath, meta, analysisMethodClass):
        self.dataPath = dataPath
        self.simName = meta.simName
        self.resDf = pd.DataFrame(np.empty((meta.shape[0],3), dtype = object), columns = ['pvalue','observed','res'], index = meta.index)
        self.metaDf = meta
        self.analysisMethodClass = analysisMethodClass
        self.analysisName = analysisMethodClass(None).methodName
        self.checkDirs()
        self.save()

    def runAnalyses(self, nperms, params = {}, distFilter = None, ba = None, clusterClient = None):
        """Run an analysis method on all the simulations, storing the results in resDf"""
        if not 'thinPermutations' in params:
            thinPermutations = 1000
        else:
            thinPermutations = params['thinPermutations']
        if not 'ditchDist' in params:
            ditchDist = True
        else:
            ditchDist = params['ditchDist']

        """Print progress to a file and to the screen"""
        with open(self.dataPath + 'sievesim/%s/%s.log' % (self.simName, self.analysisName), 'w') as logFh:
            msg = 'Running %d combinations of %d parameters' % (self.resDf.shape[0], self.resDf.shape[1] - 4)
            print msg,
            logFh.write('Starting: %s\n' % time.asctime())
            logFh.write(msg + '\n')

            startTime = time.time()
            for simID in self.resDf.index:
                """Load one simulation for analysis at a time (this clears out old ones as well)"""
                self.subLoad(simID, self.analysisName)

                tmp = analysisMethodClass(sievedata = self.metaDf.sim[simID].data)
                tmp.initialize(ba = ba, params = params)
                tmp.computeDistance()
                tmp.computeObserved(distFilter = distFilter)
                tmp.permutationTest(nperms, clusterClient)
                tmp.computePvalues()
                self.resDf.pvalue.loc[simID] = tmp.results.pvalue
                self.resDf.observed.loc[simID] = tmp.results.observed.sum()
                self.resDf.res.loc[simID] = tmp
                
                """Immediately after getting the result, save it to disk and replace with links"""
                self.subSave(simID, thinPermutations = thinPermutations, ditchDist = ditchDist)
                print '.',
                if (simID % 100) == 0:
                    logFh.write('%d: %2.1f min elapsed\n' % (simID,(time.time()-startTime)/60.))
                    logFh.flush()
            print 'done!'
            logFh.write('Finished: %s\n' % time.asctime())
    def getRelativeDir(self):
        return '%s/%s/' % (self.simName, self.analysisName)
    def getRelativePath(self, simID):
        """Used by load and save to define names of result files"""
        return self.getRelativeDir() + '%s.%d.pkl' % (self.analysisName, simID)
    def getAbsPath(self, simID):
        return self.dataPath + self.getRelativePath(simID)
    def getAbsDir(self, simID):
        return self.dataPath + self.getRelativeDir()
    def checkDirs(self):
        """Creates the needed subdirectories if they do not exist"""
        if not os.path.exists(self.getAbsDir()):
            os.makedirs(self.getAbsDir())
    def writeObject(self, obj, fn, relative = False):
        """Made this a separate function so that later I can switch from pickle easily"""
        if relative:
            absFn = self.dataPath + fn
        else:
            absFn = fn

        with open(absFn, 'wb') as fh:
            pickle.dump(obj,fh)
    def readObject(self, fn, relative = True):
        """Made this a separate function so that later I can switch from pickle easily"""
        """Fn can be relative or absolute"""
        if relative:
            absFn = self.dataPath + fn
        else:
            absFn = fn

        """Need to use the pandas function so that DataFrames are read correctly"""
        obj = pd.read_pickle(absFn)
        return obj

    def saveMetaPkl():
        self.checkDirs()

        """dict with a copy of self.resDf that contains only links to data files instead of the data"""
        columns = [c for c in self.resDf.columns if not c in ['res']]
        outDict = {'resDf':deepcopy(self.resDf[columns])}

        """Replace objects with relative filenames"""
        outDict['resDf']['res'] = outDict['resDf'].simID.map(lambda x: self.getRelativePath(x))

        """It's always safe to save this file since it should never contain data, just paths to files
           (though sometimes these data files will not exist yet at these paths (sloppy maybe?)"""
        self.writeObject(outDict, '%s/%s/%s.%s.meta.pkl' % (self.dataPath, self.simName, self.simName, self.analysisName))
    
    def _saveOne(self, simID, thinPermutations, ditchDist, preservObject = True):
        if not type(self.resDf.res[simID]) is str:
            if preserveObject:
                tmp = deepcopy(self.resDf.res.loc[simID])
            else:
                """Replace data attrib in analysisObj with filename too"""
                self.resDf.res.loc[simID].data = self.metaDf.getRelativePath(simID)
                
                """Set these attributes to None if they exist. Optionally ditch filteredDist and scannedDist too"""
                remove = ['btBA','insertBA','distFilter','dist','temp']

                if ditchDist:
                    remove += ['filteredDist','scannedDist']
                else:
                    """If don't ditch dist then keep only one of filteredDist or scannedDist depending on which exists"""
                    if hasattr(self.resDf.res[simID].results, 'scannedDist'):
                        remove += ['filteredDist']
                
                for o in remove:
                    try:
                        setattr(self.resDf.res.loc[simID].results, o, None)
                    except AttributeError:
                        pass

                if not thinPermutations is None:
                    self.resDf.res.loc[simID].thinPermutations(thinPermutations)

                tmp = self.resDf.res.loc[simID]
            
            """Save the results to a file"""
            self.writeObject(tmp, self.getRelativePath(simID))
            
            if not preserveObject:
                """Assign the object a filename link"""
                self.resDf.res.loc[simID] = self.getRelativePath(simID)

    def subSave(self, simID, thinPermutations = 1000, ditchDist = True):
        """Saves the analysis file underlying a single simulation simID to save space
           The operation replaces the object with a string to the file containing the object
           Both this function and the save function checks to make sure you are not trying to save
           a results object that is already merely a link."""
        self.saveMetaPkl()
        self._saveOne(simID, thinPermutations, ditchDist)

    def save(self, thinPermutations = 1000, ditchDist = True):
        """
        Saves the basedata, baseparams and the resDf table of all simulations to a pickled file
        The object is preserved in this process. Use subSave() when the meta set is too large to save all sims at once.

        If thinPermutations is not None then it makes a copy of resDf, thins the permutations and saves the copy
        If ditchDist = True (default) then dist and filteredDist will also be removed to save space in the file
        (leaving the existing resDf intact. if memory becomes a problem it could be the need for a deepcopy
            most of the time if we are ditching things for the file then we can probably ditch for the obj too?)
        """
        self.saveMetaPkl()

        """Go through each analysis and save the file, preserving the object in memory though"""
        for simID in self.resDf.index:
            self._saveOne(simID, thinPermutations, ditchDist)

    def loadMetaPkl(self, analysisName):
        self.resDf = self.readObject('%s/%s.%s.meta.pkl' % (self.simName, self.simName, analysisName), relative = True)['resDf']
        self.analysisName = analysisName
    def _loadOne(self, sid, deferResults, deferData):
        if not deferResults and type(self.resDf.res[sid]) is str:
            fullPath = self.dataPath + self.resDf.res.loc[sid]
            try:
            
                self.resDf.res.loc[sid] = pd.read_pickle(fullPath)
            except IOError:
                raise IOError('Result file does not exist: %s' % fullPath)

        if not deferData and type(self.resDf.res[sid].data) is str:
            fullPath = self.dataPath + self.resDf.res.loc[sid].data
            try:
                """Load the data with the relative path"""
                self.resDf.res.loc[ind].data = pd.read_pickle(fullPath)
            except IOError:
                raise IOError('Data file does not exist: %s' % fullPath)

    def load(self, analysisName, deferResults = True, deferData = True):
        """Load pickled simulation meta object and then the sieveData and sieveResults objects from the paths in resDf.sim and resDf.res
        If deferData or deferResults then only load the meta object as self.metaDf but leave paths for data/results"""

        """Load a pickled set of results from the simName and analysisName"""
        self.loadMetaPkl(analysisName)

        """Use the paths res columns to load the results into the df"""
        for ind in self.resDf.index:
            self._loadOne(ind, deferResults, deferData)

    def subLoad(self, simIDs, analysisName, deferResults = True, deferData = True):
        """Loads in the data and the results (if they exist) for a set of simulations,
           replacing the link with the actual object in resDf.
           This also replaces the old resDf object, freeing space taken by other simulations."""
        if type(simIDs) is int:
            simIDs = [simIDs]
        
        self.loadMetaPkl(analysisName)

        for sid in simIDs:
            self._loadOne(sid, deferResults, deferData)

    def eval(self, simID, evalClass = None):
        """Returns an eval object for a particular analysis
        TODO: should this file make sure that the analysis and data is loaded (not a string)?"""
        if evalClass is None:
            evalClass = siteEval
        return evalClass(self.resDf.res.loc[simID].data, self.resDf.res.loc[simID].results)

class simulationMeta(object):
    """Class for organizing, saving/loading a large number of simulations"""
    def __init__(self, dataPath, simName, basedata = None, baseparams = None, metadata = None):
        self.dataPath = dataPath
        self.simName = simName
        self.basedata = basedata
        self.baseparams = baseparams
        if not metadata is None:
            self.metaDf = metadata
    def runSimulations(self, varyParams, ba = None, nreps = 1):
        """Run a set of simulations varying the parameters over the specified ranges and storing the results in self.metaDf"""
        res = {'simID':[],'sim':[],'params':[],'rep':[]}
        res.update({vp:[] for vp in varyParams.keys()})
        
        """Add a parameter 'rep' for repeating the same simulation multiple times"""
        varyParams.update({'rep':np.arange(nreps)})
        vpKeys = varyParams.keys()

        """Initialize self.metaDf with parameters and empty sieve simulation objects"""
        for i,values in enumerate(itertools.product([varyParams[k] for k in vpKeys])):
            res['simID'].append(i)
            params = deepcopy(self.baseparams)
            for ki,k in enumerate(vpKeys):
                params.update({k:values[ki]})
                res[k].append(values[ki])
            res['params'].append(params)
            s = sieveSimulation()
            res['sim'].append(s)
        self.metaDf = pd.DataFrame(res)

        """Simulate and save as you go to conserve memory"""
        for sid in self.metaDf.index:
            if sid % round(self.metaDf.shape[0]/20) == 0:
                print "Simulation %d of %d" % (sid, self.metaDf.shape[0])
                sys.stdout.flush()
            self.metaDf.sim.loc[sid].simulate(self.basedata, params = self.metaDf.params.loc[sid], ba = ba)
            self.subSave(sid)
    def getRelativeDir(self):
        return '%s/data/' % (self.simName)
    def getRelativePath(self, simID):
        """Used by load and save to define names of result files"""
        return self.getRelativeDir() + 'simdata.%d.pkl' % (simID)
    def getAbsPath(self, simID):
        return self.dataPath + self.getRelativePath(simID)
    def getAbsDir(self, simID):
        return self.dataPath + self.getRelativeDir()
    def checkDirs(self):
        """Creates the needed subdirectories if they do not exist"""
        if not os.path.exists(self.getAbsDir()):
            os.makedirs(self.getAbsDir())
    def writeObject(self, obj, fn, relative = False):
        """Made this a separate function so that later I can switch from pickle easily"""
        if relative:
            absFn = self.dataPath + fn
        else:
            absFn = fn

        with open(absFn, 'wb') as fh:
            pickle.dump(obj,fh)
    def readObject(self, fn, relative = True):
        """Made this a separate function so that later I can switch from pickle easily"""
        """Fn can be relative or absolute"""
        if relative:
            absFn = self.dataPath + fn
        else:
            absFn = fn

        """Need to use the pandas function so that DataFrames are read correctly"""
        obj = pd.read_pickle(absFn)
        return obj

    def saveMetaPkl():
        self.checkDirs()

        """dict with a copy of self.resDf that contains only links to data files instead of the data"""
        outDict = {'basedata': self.basedata,
                   'baseparams': self.baseparams}

        """Preserve the obj so make a deepcopy of the minimal data"""
        columns = [c for c in self.metaDf.columns if not c in ['sim']]
        outDict.update({'metaDf':deepcopy(self.metaDf[columns])})

        """Replace objects with relative filenames"""
        outDict['metaDf']['sim'] = outDict['metaDf'].simID.map(lambda x: self.getRelativePath(x))

        """It's always safe to save this file since it should never contain data, just paths to files
           (though sometimes these data files will not exist yet at these paths (sloppy maybe?)"""
        self.writeObject(outDict, '%s/%s/%s.sim.meta.pkl' % (self.dataPath, self.simName, self.simName))

    def _saveOne(self, simID, preserveObject = True):
        """Save the data to a file"""
        self.writeObject(self.metaDf.sim.loc[simID], self.getRelativePath(simID))
        if not preserveObject:
            """Assign the object a filename link"""
            self.metaDf.sim[simID] = self.getRelativePath(simID)

    def subSave(self, simID, preserveObject = False):
        """Saves the analysis file underlying a single simulation simID to save space
           The operation replaces the object with a string to the file containing the object
           Both this function and the save function checks to make sure you are not trying to save
           a results object that is already merely a link."""
        self.saveMetaPkl()
        self._saveOne(simID, preserveObject)

    def save(self):
        self.saveMetaPkl()

        """Go through each data file and save it, preserving the object in memory though"""
        for simID in self.resDf.index:
            self._saveOne(simID, preserveObject = True)
    def loadMetaPkl(self, analysisName):
        tmp = self.readObject('sievesim/%s/%s.%s.meta.pkl' % (self.simName, self.simName, analysisName), relative = True)
        self.basedata = tmp['basedata']
        self.baseparams = tmp['baseparams']
        self.metaDf = tmp['metaDf']

    def _loadOne(self, sid, deferData):
        if not deferData and type(self.metaDf.sim[sid].data) is str:
            fullPath = self.dataPath + self.metaDf.sim.loc[sid].data
            try:
                """Load the data with the relative path"""
                self.metaDf.sim.loc[ind].data = pd.read_pickle(fullPath)
            except IOError:
                raise IOError('Data file does not exist: %s' % fullPath)

    def load(self, deferData = False):
        """Load pickled simulation meta object and then the sieveData and sieveResults objects from the paths in resDf.sim and resDf.res
        If deferData or deferResults then only load the meta object as self.metaDf but leave paths for data/results"""

        """Load a pickled set of results from the simName and analysisName"""
        self.loadMetaPkl()

        """Use the paths res columns to load the results into the df"""
        for ind in self.metaDf.index:
            self._loadOne(ind, deferData = deferData)

    def subLoad(self, simIDs, deferData = False):
        """Loads in the data and the results (if they exist) for a set of simulations,
           replacing the link with the actual object in resDf.
           This also replaces the old resDf object, freeing space taken by other simulations."""
        if type(simIDs) is int:
            simIDs = [simIDs]
        
        self.loadMetaPkl()

        for sid in simIDs:
            self._loadOne(sid, deferData = deferData)  

    def to_mers(self, fn = None, nmers = [9]):
        """Output all mers in all sequences in all simulations to a mers file for HLA predictions"""
        mers = set()
        for sid in self.metaDf.index:
            if type(self.metaDf.sim.loc[sid]) is str:
                self.subLoad(sid,self.analysisName)
            s = self.metaDf.sim.loc[sid]
            mers.update(set(s.to_mers(nmers = nmers, returnList = True)))
        mers = sorted(list(mers))
        if not fn is None:
            with open(fn,'w') as fh:
                for m in mers:
                    fh.write('%s\n' % m)
        else:
            return mers
    def to_hla(self, fn = None):
        """Output all HLA alleles from every participant in every simulation to a file for HLA predictions"""
        hlas = set()
        for sid in self.metaDf.index:
            if type(self.metaDf.sim.loc[sid]) is str:
                self.subLoad(sid,self.analysisName)
            s=self.metaDf.sim.loc[sid]
            hlas.update(set(s.to_hla(returnList=True)))
        hlas = sorted(list(hlas))
        if not fn is None:
            with open(fn,'w') as fh:
                for h in hlas:
                    fh.write('%s\n' % h)
        else:
            return hlas
    def saveSimPlots(self, basePath = None):
        if basePath is None:
            basePath = self.dataPath + 'sievesim/figures/'
        for simID,rec in self.metaDf.iterrows():
            figure(1)
            rec['sim'].plotEpitopes()
            savefig(basePath + '%d_sim_epitopes.png' % simID, dpi = 200)
            
            rec['sim'].plot()
            savefig(basePath + '%d_sim.png' % simID, dpi = 200)

class siteMeta(simulationMeta):
    """
    Eval functions for meta-analysis of site-based sieve analysis
    The site vs. global meta class is only important for plotting (not analysis or storage)
    """

    def plotCanon(self, variedParam = 'mutations_per_epitope', cutoff = 0.05, dx = 'vac'):
        """Scatter plot of the number of canonical vs. non-canonical sieve effects for each param value"""
        dx = 1 if dx=='vac' else 0
        mycm = ['darkgreen','fuchsia','saddlebrown','lightseagreen','gold','royalblue','tomato','thistle','tan']
    
        clf()
        uValues = sorted(self.resDf[variedParam].unique(), key = lambda x: x[dx])

        canonicalFunc = lambda aObj: ((aObj.results.observed>0) & (aObj.results.pvalue < cutoff)).sum()
        noncanonicalFunc = lambda aObj: ((aObj.results.observed<0) & (aObj.results.pvalue < cutoff)).sum()
        mx = 0
        for colori,v in enumerate(uValues):
            curInd = self.resDf[variedParam] == v
            yvec = self.resDf.res[curInd].map(canonicalFunc)
            xvec = self.resDf.res[curInd].map(noncanonicalFunc)
            mx = max(mx,max(xvec),max(yvec))
            plot(xvec, yvec, '-o', color = mycm[colori], lw = 2, label = v[dx])
            
        xlabel('Number of significant, NON-CANONICAL sites/kmers (p < %1.3f)' % cutoff)
        ylabel('Number of significant, CANONICAL sites/kmers (p < %1.3f)' % cutoff)
        xlim((0,mx+2))
        ylim((0,mx+2))
        colorLegend(colors = mycm[:colori], labels = ['%d' % v[dx] for v in uValues], title = variedParam)
        title(self.resDf.res[0].identifierStr())

def loadGlobalPvalues(dataPath, simNames, analysisNames, variedParam = 'nMutations', trtVaried = 'vaccine'):
    """Open results file for all simulation types, analysis method names and return a
       pvalue array [type x method x param x reps]"""
    
    """Used to specify which treatment group's parameter is varying"""
    trtVaried = 1 if trtVaried == 'vaccine' else 0

    getFn = lambda sname,aname: dataPath + 'sievesim/%s/%s.%s.meta.pkl' % (sname, sname, aname)

    """Load one to dims for the array"""
    sims = simulationMeta(dataPath, simNames[0])
    sims.loadFirst(getFn(simNames[0], analysisNames[0]))

    """uValues is the unique set of vac and pla parameters in tuple form (pla,vac)"""
    uValues = sorted(sims.resDf[variedParam].unique(), key = lambda param: param[trtVaried])
    paramVec = np.array([param[trtVaried] for param in uValues])

    nTypes = len(simNames)
    nMethods = len(analysisNames)
    nParams = len(uValues)
    nReps = (sims.resDf[variedParam] == uValues[0]).sum()

    pvalues = np.zeros((nTypes,nMethods,nParams,nReps), dtype = np.float64)
    for si,sname in enumerate(simNames):
        for ai,aname in enumerate(analysisNames):
            sims=simulationMeta(dataPath, sname)
            sims.loadFirst(getFn(sname, aname))

            for parami,param in enumerate(uValues):
                pvalues[si,ai,parami,:] = sims.resDf.pvalue[sims.resDf[variedParam] == param]
    return pvalues, paramVec

"""Global constants used by these plot functions"""
pFunc = lambda p: np.log10(p)
ipFunc = lambda p: 10**(p)
ypTicks = np.array([1.0, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005])
ypTickLabels = ['1.0', '0.5', '0.25', '0.1', '0.05', '0.025', '0.01', '0.005', '0.0025', '0.001', '0.0005']
xpTicks = np.array([1.0, 0.5, 0.2, 0.05, 0.02, 0.005, 0.001])
xpTickLabels = ['1.0', '0.5', '0.2', '0.05', '0.02', '0.005', '0.001',]

def plotPRByP(sname,simNames,analysisNames,pvalues,parami,colors,labels):
    """Plot positive rate (FPR or TPR depending on sim parameters) for each method as a function of the pvalue"""
    nPerms=ceil(1./pvalues.min())
    nTypes,nMethods,nParams,nReps=pvalues.shape

    si=simNames.index(sname)
    """
    yvecs=zeros((len(analysisNames),len(pThresh)))
    for ai,aname in enumerate(analysisNames):
        for threshi,thresh in enumerate(pThresh):
            yvecs[ai,threshi]=(pFunc(pvalues[si,ai,parami,:]) > pThresh[threshi]).sum()/nReps
    """

    cla()
    plot(pFunc(np.array([0.0001,1.2])),pFunc(np.array([0.0001,1.2])),'--',color='gray')
    for ai,aname in enumerate(analysisNames):
        #plot(pThresh,pFunc(yvecs[ai,:]),'-',color=colors[aname],label=labels[aname],lw=2)
        plot(sort(pFunc(pvalues[si,ai,parami,:])),pFunc(arange(nReps)/nReps),'-',color=colors[aname],label=labels[aname],lw=2)
    xlabel('Significance threshold (p)',size='x-large')
    #ylabel('Detection rate',size='x-large')
    xt=xticks()
    #xticks(xt[0],ipFunc(xt[0]),size='large')
    xticks(pFunc(xpTicks),xpTickLabels,size='large')
    yticks(pFunc(ypTicks),ypTickLabels,size='large')
    
    xlim(0.001,1.1)
    ylim(0.001,1.1)

    #ylim(pFunc(array([1/nReps,1.0])))
    """Reverse the x-axis"""
    #xlim(pFunc(array([0.001,0.2])))
    #xlim((-1,round(pFunc(1/nPerms))+1))
    #title('%s: param = %d' % (labels[sname],parami),size='x-large')
    title('%s' % (labels[sname]),size='x-large')
    #legend(loc='upper right')

def plotFPR(analysisNames,pvalues,colors,labels):
    """Plot false-positive rate (FPR) for each method as a function of the pvalue
       Uses the no_mutations simName"""
    
    nTypes,nMethods,nParams,nReps=pvalues.shape
    nReps=nReps*nTypes
    nPerms=ceil(1./pvalues.min())

    parami=0

    clf()
    plot(pFunc(np.array([0.001,0.2])),pFunc(np.array([0.001,0.2])),'--',color='gray')
    for ai,aname in enumerate(analysisNames):
        plot(sort(pFunc(pvalues[:,ai,parami,:].flatten())),pFunc(arange(nReps)/nReps),'-',color=colors[aname],label=labels[aname],lw=2)
    xticks(pFunc(xpTicks),xpTickLabels,size='large')
    yticks(pFunc(ypTicks),ypTickLabels,size='large')
    
    ylim(pFunc(np.array([1/nReps,1.0])))
    xlim(pFunc(np.array([0.001,0.2])))
    xlabel('Significance threshold (p)',size='x-large')
    ylabel('False discovery rate (FDR)',size='x-large')


def plotPRatP(sname,simNames,analysisNames,pvalues,paramVec,paramLabel,colors,labels,cutoff=0.05):
    """Plot positve rate as a function of the paramVec given significance cutoff"""
    nPerms=ceil(1./pvalues.min())
    nTypes,nMethods,nParams,nReps=pvalues.shape

    si=simNames.index(sname)
    yvecs=((pvalues<cutoff).sum(axis=3)/nReps)[si,:,:]

    cla()
    for ai,aname in enumerate(analysisNames):
        plot(paramVec,yvecs[ai,:],'-o',color=colors[aname],label=labels[aname],lw=2)
    xlabel(paramLabel,size='medium')
    #ylabel('Detection rate',size='x-large')
    xticks(paramVec,['%1.0f' % p for p in paramVec],size='large')
    yticks(size='large')
    ylim((-0.05,1.05))
    xlim((paramVec[0]-1,paramVec[-1]+1))
    title('%s' % (labels[sname]),size='x-large')
    #legend(loc='upper right')


"""
NOTE: these functions won't work anymore now that I've changed around the loading and unloading
"""
def plotVariedParamGlobal(dataPath,analysisMethodNames,simName,variedParam='mutations_per_epitope',dx='vac'):
    """Plot the observed magnitude and significance of sieve analysis over the range of the varied parameter"""
    dx = 1 if dx=='vac' else 0
    mycm = ['darkgreen','fuchsia','saddlebrown','lightseagreen','gold','royalblue','tomato','thistle','tan']
    cutoff = 0.05

    clf()
    colori=0
    for ai,aName in enumerate(analysisMethodNames):
        if aName.find('global')>=0:
            outFn=dataPath + 'sievesim/%s/%s.%s.meta.pkl' % (simName,simName,aName)
            sims=simulationMeta()
            print 'Loading %s...' % outFn,
            sims.load(outFn)
            print 'done'

            uValues=sorted(sims.resDf[variedParam].unique(),key=lambda x: x[dx])

            xvec=[]
            yvec=[]
            for v in uValues:
                xvec.append(v[dx])
                yvec.append(sims.resDf.res[sims.resDf[variedParam]==v].map(lambda aObj: aObj.results.observed[0]).mean())
            subplot(2,1,1)
            plot(xvec,yvec,'-o',color=mycm[colori],lw=2,label=aName)

            xvec=[]
            yvec=[]
            for v in uValues:
                """
                tmpy=list(sims.resDf.res[sims.resDf[variedParam]==v].map(lambda aObj: -10*np.log10(aObj.results.pvalue[0])))
                yvec.extend(tmpy)
                xvec.extend([v[dx]]*len(tmpy))
                """
                xvec.append(v[dx])
                tmpy=sims.resDf.res[sims.resDf[variedParam]==v].map(lambda aObj: aObj.results.pvalue[0])
                #tmpy[tmpy<(1/1e4)]=1/1e4
                yvec.append(tmpy.map(pFunc).mean())

            subplot(2,1,2)
            plot(xvec,yvec,'-o',color=mycm[colori],lw=2,label=aName)
            colori+=1

    subplot(2,1,1)
    #xlabel(variedParam)
    ylabel('Observed t-statistic')
    legend(loc='best')
    title('%s\nSimulation: %s' % (sims.resDf.res[0].identifierStr(),simName))

    subplot(2,1,2)
    axhline(pFunc(cutoff), ls = '--', color = 'gray')
    #plot(xlim(),[-10*log10(cutoff)]*2,'--',color='gray')
    xlabel(variedParam)
    #ylabel('$-10 log_{10}(p)$')
    pLabels = np.array([1,0.2,0.1,0.05,1e-2,1e-3,1e-4])
    yticks(pFunc(pLabels), pLabels)
    ylim((0,45))
    ylabel('Permutation test p-value')


def plotVariedParamSite(dataPath,analysisMethodNames,simName,variedParam='mutations_per_epitope',dx='vac'):
    """Plot the observed magnitude, significance, specificity and sensitivity of sieve analysis 
    over the range of the varied parameter"""
    dx = 1 if dx=='vac' else 0
    cutoff = 0.05
    mycm = ['darkgreen','fuchsia','saddlebrown','lightseagreen','gold','royalblue','tomato','thistle','tan']
    
    clf()
    colori = 0
    for ai,aName in enumerate(analysisMethodNames):
        if aName.find('site')>=0 or aName.find('kmer')>=0:
            outFn=dataPath + 'sievesim/%s/%s.%s.meta.pkl' % (simName,simName,aName)
            sims=simulationMeta()
            sims.load(outFn)

            uValues=sorted(sims.resDf[variedParam].unique(),key=lambda x: x[dx])
            xvec=[v[dx] for v in uValues]
            yvec=[]
            for v in uValues:
                yvec.append(sims.resDf.res[sims.resDf[variedParam]==v].map(lambda aObj: (aObj.results.pvalue<cutoff).sum()).mean())
            plot(xvec,yvec,'-o',color=mycm[colori],lw=2,label=aName)
            colori+=1
            
    ylabel('Number of significant sites/kmers (of %d otal)' % (len(sims.resDf.res[0].data.insertSeq)))
    legend(loc='best')
    title('%s\nSimulation: %s' % (sims.resDf.res[0].identifierStr(),simName))
    xlabel(variedParam)
    xticks(xvec)
    xlim((min(xvec)-1,max(xvec)+1))

def plotSiteROC(sims,variedParam='mutations_per_epitope',cutoff=0.05,dx='vac'):
    """Scatter plot of the specificity and sensitivity of detecting the simulated Ab epitope (mutation) sites
       (only makes sense for a vxmatch_siteAnalysis of Ab epitopes)
       [('A_6801', 'ATLYCVHQK', (82, 83, 84, 85, 86, 87, 88, 89, 90), 5.849),...]
       [('Ab', 'G', array([23]), nan),...]
       []
    """
    dx = 1 if dx=='vac' else 0
    mycm=['darkgreen','fuchsia','saddlebrown','lightseagreen','gold','royalblue','tomato','thistle','tan']

    clf()
    uValues = sorted(sims.resDf[variedParam].unique(), key = lambda x: x[dx])
    yvec = np.empty(len(uValues),dtype = np.float64)
    xvec = np.empty(len(uValues),dtype = np.float64)
    for colori,v in enumerate(uValues):
        """inds for the reps matching param v"""
        curInd = sims.resDf[variedParam]==v
        tprArr=empty(curInd.sum(),dtype=float64)
        fprArr=empty(curInd.sum(),dtype=float64)
        for indi,ind in enumerate(where(curInd)[0]):
            aObj = sims.resDf.res[ind]
            allSites = set(arange(len(aObj.data.insertSeq)))
            abSites = set([e[2][0] for e in aObj.data.simDf.epitopes[0]])
            nonAbSites = allSites.difference(abSites)
            
            posSites = set(where((aObj.results.observed>0) & (aObj.results.pvalue<cutoff))[0])
            negSites = set(where((aObj.results.observed<0) | (aObj.results.pvalue>cutoff))[0])
            
            truePosSites = posSites.intersection(abSites)
            trueNegSites = negSites.intersection(nonAbSites)
            falsePosSites = posSites.intersection(nonAbSites)
            
            tprArr[indi] = len(truePosSites)/len(abSites)
            fprArr[indi] = len(falsePosSites)/len(allSites)
        yvec[colori] = tprArr.mean()
        xvec[colori] = fprArr.mean()
    plot(xvec, yvec, '-o', color = mycm[0], lw = 2)
        
    xlabel('False positive rate')
    ylabel('True positive rate')
    xlim((0,1))
    ylim((0,1))
    title('%s\nDetecting simulated Ab contact sites (p < %1.3f)' % (sims.resDf.res[0].identifierStr(),cutoff))
