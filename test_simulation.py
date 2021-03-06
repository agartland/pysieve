import unittest
import numpy as np

from simulation import *
from za_hla import loadSouthAfricaHLA
from HLAPredCache import hlaPredCache, RandCache
from meta import *

from analysis_substbased import *

from seqdistance.matrices import addGapScores, binarySubst, binGapScores

import ipdb

class TestSimulation(unittest.TestCase):
    def setUp(self):
        pass
    def getBaseSimulation(self):
        freq = loadSouthAfricaHLA('./data/LANL_study4_patients.csv')
        simSD = simSDFromLANL('./data/', protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
        s = sieveSimulation()
        return s, simSD
    def test_parseHLA(self):
        freq = loadSouthAfricaHLA('./data/LANL_study4_patients.csv')
    
    def test_loadingSeedData(self):
        freq = loadSouthAfricaHLA('./data/LANL_study4_patients.csv')
        simSD = simSDFromLANL('./data/', protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
        self.assertTrue(hasattr(simSD,'hlaFreq'))
    
    def test_oneInsertBased(self):
        params =  { 'n':(10, 10),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'escapeMutations':(False,False),
                    'insertAsBase':(True,True)}       #Use vaccine insert

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params)
        self.assertEqual(s.data.insertSeq, simSD.insertSeq)
        self.assertEqual(s.data.seqDf.shape[0], params['n'][0] + params['n'][1])
        self.assertTrue(s.data.simParams == params)
        #ipdb.set_trace()
        self.assertEqual(np.median(s.data.simDf['epitopes'].map(len)[s.data.ptidDf.vaccinated]), params['nEpitopes'][1] * 9)
        #self.assertEqual(np.median(s.data.simDf['mutations'].map(len)[s.data.ptidDf.vaccinated]), params['nMutations'][1])
    
    def test_oneAb(self):
        params =  { 'n':(10, 10),
                    'nMutations': (0, 4), 
                    'nEpitopes': (0, 3),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params)
        self.assertEqual(s.data.insertSeq, simSD.insertSeq)
        self.assertEqual(s.data.seqDf.shape[0], params['n'][0] + params['n'][1])
        self.assertTrue(s.data.simParams == params)
        
        self.assertEqual(np.median(s.data.simDf['epitopes'].map(len)[s.data.ptidDf.vaccinated]), params['nEpitopes'][1] * 9)
        #self.assertEqual(np.median(s.data.simDf['mutations'].map(len)[s.data.ptidDf.vaccinated]), params['nMutations'][1])
    
    def test_oneNoEpitopes(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 0),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params)

        self.assertEqual(s.data.insertSeq, simSD.insertSeq)
        self.assertEqual(s.data.seqDf.shape[0], params['n'][0] + params['n'][1])
        self.assertTrue(s.data.simParams == params)

    def test_meta(self):
        s,simSD = self.getBaseSimulation()
        varyParams = {'nMutations' : [(0,n) for n in [3,6,9]]}

        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        sims = simulationMeta('./data/', 'test_simulations', basedata = simSD, baseparams = params, verbose = False)
        sims.runSimulations(varyParams, ba = None, nreps = 5)
        self.assertEqual(sims.metaDf.shape[0], 5 * 3)
        self.assertTrue(isinstance(sims.metaDf.sim.iloc[0], basestring))
        sims.save()
        self.assertTrue(isinstance(sims.metaDf.sim.iloc[0], basestring))
        sims.subLoad([0])
        print type(sims.metaDf.sim.iloc[0])
        self.assertTrue(isinstance(sims.metaDf.sim.iloc[0], sieveSimulation))
        
    @unittest.skip('d')
    def test_meta_analysis(self):
        s,simSD = self.getBaseSimulation()
        varyParams = {'nMutations' : [(0,n) for n in [3,6,9]]}

        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        sims = simulationMeta('./data/', 'test_simulations', basedata = simSD, baseparams = params, verbose = False)
        sims.runSimulations(varyParams, ba = None, nreps = 2)
        ma = analysisMeta('./data/', sims, vxmatch_globalAnalysis, verbose = False)
        analysisParams = {'nmer':9,
                          'subst': addGapScores(binarySubst, binGapScores)}
        ma.runAnalyses(50, params = {}, distFilter = None, ba = None, clusterClient = None)
        self.assertEqual(ma.resDf.shape[0], sims.metaDf.shape[0])

class TestTcellEpitopeSimulations(unittest.TestCase):
    def setUp(self):
        fn = './data/gag.cladec.ZA.2012.netmhcpan'
        #self.ba = hlaPredCache(fn, kmers = [9])
        self.ba = RandCache()
    def getBaseSimulation(self):
        freq = loadSouthAfricaHLA('./data/LANL_study4_patients.csv')
        simSD = simSDFromLANL('./data/', protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
        s = sieveSimulation(testMode = True)
        return s,simSD
    def test_noBA(self):
        params =  { 'n':(15, 15),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,False),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        with self.assertRaises(ValueError) as cm:
            s.simulate(simSD, params = params, ba = None)

    def test_oneTcellEpitopes(self):
        params =  { 'n':(15, 15),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,False),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params, ba = self.ba)
        self.assertEqual(s.data.insertSeq, simSD.insertSeq)
        self.assertEqual(s.data.seqDf.shape[0], params['n'][0] + params['n'][1])
        self.assertTrue(s.data.simParams == params)
        #ipdb.set_trace()
        self.assertEqual(np.median(s.data.simDf['epitopes'].map(len)[s.data.ptidDf.vaccinated]), params['nEpitopes'][1])
        #self.assertEqual(np.median(s.data.simDf['mutations'].map(len)[s.data.ptidDf.vaccinated]), params['nMutations'][1])

    def test_oneTcellEscape(self):
        params =  { 'n':(15, 15),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,True),
                    'escapeDelta':(0,2),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params, ba = self.ba)
        self.assertEqual(s.data.insertSeq, simSD.insertSeq)
        self.assertEqual(s.data.seqDf.shape[0], params['n'][0] + params['n'][1])
        self.assertTrue(s.data.simParams == params)
        #ipdb.set_trace()
        self.assertEqual(np.median(s.data.simDf['epitopes'].map(len)[s.data.ptidDf.vaccinated]), params['nEpitopes'][1])
        #self.assertEqual(np.median(s.data.simDf['mutations'].map(len)[s.data.ptidDf.vaccinated]), params['nMutations'][1])
    
    @unittest.skip("Plotting not implemented")
    def test_plotSimulation(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,False),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s,simSD = self.getBaseSimulation()
        s.simulate(simSD, params = params)

        s.plot()
        s.plotEpitopes()

if __name__ == '__main__':
    unittest.main()