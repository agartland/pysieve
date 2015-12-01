import numpy as np
import pandas as pd

from HLAPredCache import hlaPredCache, RandCache
import pysieve
from pysieve import tcell, substbased
from vtn_sieve import sieveVTN503, sieveVTN502, sieveVTN505, sieveRV144
from seqdistance.matrices import binarySubst, addGapScores, binGapScores

class sieveHIV(pysieve.sieveDataMethods):
    _validProteinInsertNames={'gag':'MRK_INSERT_GAG',
                              'nef':'MRK_INSERT_NEF',
                              'pol':'MRK_INSERT_POL'}

    @property
    def validAnalyses(self):
        return [(k,self._validProteinInsertNames[k]) for k in ['gag','pol','nef']]

    def loadData(self, dataPath=None, proteinName='gag', insertName='NA', regionInds=None):
        """Loads breakthrough sequences, vaccine sequences, HLA alleles and other participant data from study-specific files."""
        pass


"""Instantiation makes available some basic info about the trial
(e.g. a list of protein and insert names that are valid)"""
base = sieveVTN503()
validAnalyses = base.validAnalyses
print base.studyName

analysisParams = {'subst':addGapScores(binarySubst, binGapScores)}

for va in base.validAnalyses:
    print va['proteinName'], va['insertName']

s = sieveVTN503()
s.loadData(proteinName='gag', insertName='MRK')

a = pysieve.analysis_substbased.vxmatch_globalAnalysis(sievedata=s.data)

a.initialize(params=analysisParams)
a.computeDistance(params=analysisParams)
site_filter = pysieve.filters.diversityFilter(s.data, minMM=5, minM=5)
a.computeObserved(distFilter=site_filter)

a.permutationTest(10000)
a.computePvalues()

resDf = a.to_df()

"""e = globalEval(a.data, a.results)
figh = plt.figure(1)
e.plotUnblindedDistance(figH)
figh.savefig(DATA_PATH + 'figures/%s.%s.unblinded_dist.png' % (proteinName,insertName))"""