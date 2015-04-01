from pylab import *
from pysieve import *

class sieveHIV(sieveDataMethods):
    _validProteinInsertNames={'gag':'MRK_INSERT_GAG',
                              'nef':'MRK_INSERT_NEF',
                              'pol':'MRK_INSERT_POL'}

    @property
    def validAnalyses(self):
        return [(k,self._validProteinInsertNames[k]) for k in ['gag','pol','nef']]

    def loadData(self,dataPath=None,proteinName='gag',insertName='NA',regionInds=None):
        """Loads breakthrough sequences, HLA alleles and other participant data from files."""
        pass


base = sieveHIV()
validAnalyses = base.validAnalyses

analysisParams = {'nmer':9,
                  'binding':5,
                  'escape':7,
                  'minDelta':1,
                  'ignoreGappedKmers':False,
                  'subst':addGapScores(binarySubst,binGapScores)}

figH = figure(1)

for proteinName,insertName in base.validAnalyses:
    data = sieveHIV()
    data.loadData(proteinName=proteinName, insertName=insertName)

    site_filter = diversityFilter(d, minMM=5, minM=5)

    a = vxmatch_globalAnalysis(sievedata = data)

    a.initialize(params = analysisParams)

    a.computeDistance()
    a.computeObserved(filter = site_filter)

    a.permutationTest(10000, remoteClient = rclient)
    a.computePvalues()
    
    dataFn,analysisFn = saveSieve(DATA_PATH,a)
    a.to_csv(DATA_PATH + '%s.%s.results.csv' % (proteinName,insertName))

    e = globalEval(a.data, a.results)

    e.plotUnblindedDistance(figH)
    figH.savefig(DATA_PATH + 'figures/%s.%s.unblinded_dist.png' % (proteinName,insertName))

