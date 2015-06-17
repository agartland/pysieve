import unittest
import numpy as np

from simulation import *
from za_hla import loadSouthAfricaHLA
from hla_prediction import hlaPredCache

class TestSimulation(unittest.TestCase):
    def setUp(self):
        fn = DATA_PATH + 'data/gag.cladec.ZA.2012.netmhcpan'
        self.ba = hlaPredCache(fn,kmers=[9])
        #But insert predictions will be missing. Load a full file for testing instead"""
    def getBaseSimulation(self):
        freq = loadSouthAfricaHLA()
        simSD = simSDFromLANL(protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
        s = sieveSimulation()
        return s
    def test_parseHLA(self):
        freq = loadSouthAfricaHLA()
    def test_loadingSeedData(self):
        freq = loadSouthAfricaHLA()
        simSD = simSDFromLANL(protein = 'gag', year = 2012, clade = 'C', hlaFreq = freq)
    def test_oneInsertBased(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'escapeMutations':(False,False),
                    'insertAsBase':(True,True)}       #Use vaccine insert

        s = self.getBaseSimulation()
        s.simulate(s, params = params)    

    def test_oneAb(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s = self.getBaseSimulation()
        s.simulate(s, params = params)
    def test_oneNoEpitopes(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 0),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s = self.getBaseSimulation()
        s.simulate(s, params = params)
    def test_oneTcellEpitopes(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,False),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s = self.getBaseSimulation()
        s.simulate(s, params = params)
    def test_oneTcellEscape(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,True),
                    'escapeDelta':(0,2),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s = self.getBaseSimulation()
        s.simulate(s, params = params)
    def test_meta(self):
        s = self.getBaseSimulation()
        varyParams = {'nMutations' : [(0,n) for n in [3,6,9]]}

        params =  { 'n':(5, 5),
                    'nMutations': (0, 3), 
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (None,None),      #None means Ab epitope for pickEpitopes()
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        sims = simulationMeta('data/', 'test_simulations', basedata = s, baseparams = params)
        sims.runSimulations(varyParams, ba = self.ba, nreps = 5)
        sims.save()
    def test_loadMeta(self):
        pass
    def test_plotSimulation(self):
        params =  { 'n':(5, 5),
                    'nMutations': (0, 3),
                    'nEpitopes': (0, 2),
                    'epitopeThreshold': (0,5),
                    'escapeMutations':(False,False),
                    'insertAsBase':(False,False)}       #Use sampled sequences as base

        s = self.getBaseSimulation()
        s.simulate(s, params = params)

        s.plot()
        s.plotEpitopes()



