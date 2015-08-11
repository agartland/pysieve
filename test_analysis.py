import unittest
import numpy as np

from simulation import *
from za_hla import loadSouthAfricaHLA
from HLAPredCache import hlaPredCache
from meta import *
from analysis_substbased import *
from analysis_tcell import *

from seqdistance.matrices import addGapScores, binarySubst, binGapScores
import ipdb

class TestSubstAnalysis(unittest.TestCase):
    def setUp(self):
        params =  { 'n':(40, 40),
                    'nMutations': (0, 4), 
                    'nEpitopes': (0, 3),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params, verbose = False)
        self.sim = s

    def getBaseSimulation(self):
        freq = loadSouthAfricaHLA('./data/LANL_study4_patients.csv')
        simSD = simSDFromLANL('./data/', protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
        s = sieveSimulation()
        return s, simSD
    
    def test_vxmatch_global(self):
        analysisParams = {'nmer':9,
                          'subst': addGapScores(binarySubst, binGapScores)}
        a = vxmatch_globalAnalysis(sievedata = self.sim.data)

        a.computeDistance(params = analysisParams)
        
        a.computeObserved()

        a.permutationTest(100, clusterClient = None)
        a.computePvalues()
        
        #dataFn, analysisFn = saveSieve('./data/',a)
        a.to_csv('./data/test_vxmatch.results.csv')

if __name__ == '__main__':
    unittest.main()