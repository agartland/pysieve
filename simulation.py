from __future__ import division
"""
simulation.py
contains classes for creating sieve datasets

TODO: it would be nice to have a mode that generated Ab-based pressure without needing to load HLA data in anyway.
I want a oneliner to generate sieve datasets for experimenting with ML methods...

common input is a sieveData object with:
    (1) insert sequence
    (2) base sequences (use placebo)
    (3) parameters
    (4) HLA distribution (other ptid covariates?)
common methods:
    (1) read inputs from files
    (2) simulate variation, given parameters/inputs
    (3) store sieve dataset in a class
    (4) store sieve dataset in a FASTA or CSV file

Algorithm:

For each person in the study:

Start with sequence sampled from the seqDf (sampled w/o replacement)
Assign HLA alleles based on frequencies specified
Choose [N_eptitopes] based on [epitope_threshold] and a random HLA allele of the participant

Repeat [N_eptitopes] x [mutations_per_epitope] times:

    * Choose a site at random with probability [probability_mutation_in_epitope] that the site is in a pre-specified epitope
    * Mutate site by choosing an AA from possible AA defined by the supplied seqDf. If [mutation_rejection] then reject mutations until one is found that decreases current binding affinity of that epitope


Method analysis:
Test out global, k-mer specific and site-specific methods on simulated data sets
Test power to reject null hypothesis (H0: V=P) under various parameter combinations:

    * Systematically increasing the sieve effect by increasing mutations_per_epitope
    * Varying breadth (N_epitopes)
    * Varying number of non-epitope related mutations
    * Alter whether or not all mutations reduce binding affinity

One advantage of mutating sites rather than k-mers is that there is a potential for epitopes to overlap and compete/cancel with eachother.

Would it be better just to simulate shifts in binding affinity for a certain number of epitopes? Placebo recipients would then just have either: 

    1. a higher epitope threshold
    2. fewer mutations (affinity shifts) per epitope
    3. or same epitopes but neutral affinity shifts
"""
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import pandas as pd
#from Queue import Queue,Empty
#import threading
from numpy import *
from random import choice
from numpy.random import permutation, randint
from copy import deepcopy
from scipy import stats
from seqplot import *
from hla_prediction import *
from discreteSampler import discreteSampler

try:
    import mhc_bindings.mhc_bindings as mb
except ImportError:
    print 'Could not load mhc_bindings for sequence simulation.'
    
from data import *
from sieveio import *
from za_hla import loadSouthAfricaHLA

__all__ = ['sieveSimulationData',
          'sieveSimulation',
          'simSDFromLANL']

class sieveSimulationData(sieveData):
    """Expanded sieveData class to hold simulation metadata"""
    hlaAFreq=None
    hlaBFreq=None
    simParams=None
    date=None
    """Columns for each sequence: mutations, insertEpitopes"""
    simDf=None

class sieveSimulation(sieveDataMethods):
    data=None
    
    def __init__(self,sievedata=sieveSimulationData(),predMethod='netMHCpan',testMode=False):
        self.data=deepcopy(sievedata)
        """The Q for doing predictions during the simulation"""
        #self.requestQ=Queue()
        self.predMethod=predMethod
        self.testMode=testMode

    def startFetcher(self):
        """Will only need to start the fetcher if its in escape mode"""
        self.fetcherThread=threading.Thread(target=self.fetcher)
        self.fetcherThread.daemon=True
        self.fetcherThread.start()

    def stopFetcher(self):
        self.requestQ.put(('STOP','STOP','STOP'))

    def fetcher(self):
        """Centralized predictor thread.
        From reqQ, gets an HLA allele, a list of peptide variants and an outgoing Q
        Run all predictions and return the results in a pd.DataFrame() with columns hla, peptide, ic50 and index of ints"""
        count=1
        while True:
            """Effectively a blocking call to get next request"""
            try:
                hla,peptides,outQ=self.requestQ.get(block=False,timeout=2)
                if hla=='STOP':
                    print 'Fetcher is dying...'
                    break

                """Do the predictions (they might be slow)"""
                df=self.predictHLABinding(hla,peptides)
                outQ.put(df)
                print 'Fetcher made decision number %d' % count
                count+=1
            except Empty:
                pass
    def predictHLABinding(self,hla,peptides):
        """Run all predictions and return the results in a pd.DataFrame()
        with columns hla, peptide, ic50 and index of ints
        Predictions will not be in the same order as requested!
        Prediction is I/O bound so may be appropriate to put it in a fetcher thread?"""
        if not self.testMode:
            """See if its in the self.ba. If not then get the prediction and add it to self.ba. Pull all results from self.ba"""
            neededPeptides=[p for p in peptides if not self.ba.has_key((hla.replace('*','_'),p))]
            if len(neededPeptides)>0:
                results=mb.fetch_return(methods=[self.predMethod],
                                        mhc_strs=[hla],
                                        peptide_strs=neededPeptides,
                                        peptide_length=9,
                                        force=True,
                                        query=False,
                                        verb=False,
                                        no_cache=True,
                                        db_passwd='R1ghtN0w',login='agartlan')
                addHLAs,addPeptides,addValues=[],[],[]
                for method,hla,pep,core,pred in results:
                    addPeptides.append(pep)
                    addHLAs.append(hla)
                    addValues.append(pred)
                self.ba.addPredictionValues(addHLAs,addPeptides,addValues)
            results=[(self.predMethod,hla,p,'CORE',self.ba[(hla,p)]) for p in peptides]
        else:
            results=[('DUMMY',hla,pep,pep,rand()*15) for p in peptides]

        d={'hla':[],'peptide':[],'ic50':[]}
        for method,hla,pep,core,pred in results:
            d['hla'].append(hla)
            d['peptide'].append(pep)
            d['ic50'].append(pred)
        df=pd.DataFrame(d)
        return df

    def simulate(self,sd,params=None,ba=None):
        """Simulate a sieve dataset based on the sieve data object sd and given params.
        Within the sd the following are required:
            insertSeq
            hlaDf (or hlaFreq specified)
            ptidDf - will limit to placebo sequences (column: ~vaccinated)
            HXB2
            mapDf
            studyName
            proteinName
        """

        """All params have a vaccine (1) and placebo (0) value (placebo, vaccine)"""
        if params is None:
            raise Exception('No params[] specified!')
        if ba is None:
            """Load binding affinity dict for identifying insert epitopes"""
            ba=hlaPredCache(params['baFn'],kmers=[9])
        self.ba=ba
        self.ba.warn=True
        """If we'll need to use the predictor, get it started (for now no, trying prediction without multithreading"""
        """
        if any(params['escapeMutations']):
            self.startFetcher()
        """

        """Take the HLA frequencies from the participants or from specified frequencies in sd.hlaFreq['A']"""
        if not hasattr(sd,'hlaFreq'):
            """Limit HLA distribution to uninfected and placebo infected participants"""
            hlaDf=sd.hlaDf.join(sd.ptidDf)
            hlaDf=hlaDf[(~hlaDf.infected) | (~hlaDf.vaccinated & hlaDf.infected)]

            uHLA4_A=[]
            uHLA4_B=[]
            for h in sd.uHLA4:
                if h[0]=='A':
                    uHLA4_A.append(h)
                elif h[0]=='B':
                    uHLA4_B.append(h)

            hlaAFreq=hlaDf[uHLA4_A].sum()/hlaDf[uHLA4_A].sum().sum()
            hlaBFreq=hlaDf[uHLA4_B].sum()/hlaDf[uHLA4_B].sum().sum()
        else:
            hlaAFreq=sd.hlaFreq['A']
            hlaBFreq=sd.hlaFreq['B']

        def pickHLAs(hlaAFreq,hlaBFreq):
            tmpA=discreteSampler(hlaAFreq,2)
            tmpB=discreteSampler(hlaBFreq,2)
            return tmpA+tmpB
        #randHLAA=discreteSampler(hlaAFreq,2*2*(params['n'][0]+params['n'][1]))
        #randHLAB=discreteSampler(hlaBFreq,2*2*(params['n'][0]+params['n'][1]))
        #hlaRandi=0

        """Limit sequences to placebo sequences."""
        baseSeq=sd.seqDf.join(sd.ptidDf)
        baseSeq=baseSeq[~baseSeq.vaccinated]['seq']

        """Make a list of possible AA for each site."""
        baseSeqMat=align2mat(baseSeq)
        possibleAA=[''.join(unique(baseSeqMat[:,coli])).replace('-','') for coli in arange(baseSeqMat.shape[1])]
        
        """Only sites that have no seqs with a gap and that have more than 1 possible AA are valid epitope indices (ie sites for a mutation)"""
        """validEpitopeInds=permutation(find(~any(baseSeqMat=='-',axis=0) & (array([len(p)>1 for p in possibleAA]))))"""

        """All sites that do not have a consensus gap and that have more than 1 possible AA, are possible mutations sites"""
        validEpitopeInds=permutation(find(array([aa!='-' for aa in sd.insertSeq]) & array([len(p)>1 for p in possibleAA])))
        #print "There are %d sites possible for mutations" % len(validEpitopeInds)

        alphabet=''.join(set(sd.insertSeq.__str__()))

        dataColumns=['ptid','vaccinated','hla','seq','seqID','epitopes','mutations','infected']
        data={c:[] for c in dataColumns}

        """Permutation of seq indices so that we can pick seed sequences without replacement"""
        seqInds=permutation(baseSeq.index)
        
        ptid=0
        #print "Simulating %d vaccine and %d placebo recipients' isolates..." % (params['n'][1],params['n'][0])
        for vaccinated in [True,False]:
            for i in arange(params['n'][vaccinated]):
                """
                if vaccinated:
                    print 'V',
                else:
                    print 'P',
                """
                """Pick sequence at random from the base set"""
                if not params['insertAsBase'][vaccinated]:
                    basePtid=seqInds[ptid]
                    seq=baseSeq[basePtid]
                else:
                    seq=sd.insertSeq
                    basePtid='insert'

                """List to keep track of mutations"""
                mutations=[]
                
                epitopeInds,epitopes=[],[]
                """If the chosen HLA alleles lead to no potential epitopes that are less than the threshold..."""
                while len(epitopeInds) == 0:
                    tmpHLA=pickHLAs(hlaAFreq, hlaBFreq)
                    """Returns a tuple of indices where each set of indices identify AA locations that make up a single epitope"""
                    epitopeInds,epitopes=pickEpitopes(sd.insertSeq,seq,params['nEpitopes'][vaccinated],params['epitopeThreshold'][vaccinated],ba,tmpHLA,validEpitopeInds)
                    if len(epitopeInds) == 0:
                        print 'No epitopes for %s, %s, %s, %s' % tuple(tmpHLA)
                    else:
                        if params['escapeMutations'][vaccinated] and any([seq[e[2][0]]=='-' for e in epitopes]):
                            epitopeInds,epitopes=[],[]
                            """If the first AA of the epitope in the BT is - then find a different epitope"""
                
                '''
                """For each site in an epitope induce a mutation with Pr = epitopeMutation * (nMutations/len(epitopeInds))
                This algorithm allows for 1 or 0 mutations at each site (never 2 mutations at 1 site)
                If epitopeMutation is less than 1 it does not mean mutations will be introduced elsewhere
                (theoretically the probability could be spread evenly across non epitope sites, but this is for another day)"""
                if len(epitopeInds)<1:
                    probMut=0.
                else:
                    probMut=params['epitopeMutation'][vaccinated]*params['nMutations'][vaccinated]/float32(len(epitopeInds))
                for sitei in epitopeInds:
                    if rand()<probMut:
                        """Mutate site based on circulating or breakthrough strains or randomly from alphabet"""
                        #newAA=choice(alphabet)
                        tmpPossibleAA=deepcopy(possibleAA[sitei]).replace(seq[sitei],'')
                        if len(tmpPossibleAA)==0:
                            raise Exception('No possible AA for site %d' % sitei)
                        newAA=choice(tmpPossibleAA)
                        mutations.append((sitei,seq[sitei],newAA))
                        seq=seq[:sitei] + newAA + seq[sitei+1:]
                        """ADD SOMETHING HERE ABOUT FORCING ESCAPE MUTATIONS"""
                '''
                """NEW ALGORITHM
                Use binomial dist to first determine how many mutations are to be made
                while there are still mutations to be made:
                    Pick an epitope at random
                    Pick a site from the epitope at random (or pick epitope variants of 1 mut)
                    Mutate the site (or accept the escape variant)
                    If was acceptable, increment the counter"""
                if len(epitopeInds)>=1:
                    probMut=params['nMutations'][vaccinated]/float32(len(epitopeInds))
                    randVar=stats.binom(n=len(epitopeInds),p=probMut)
                    actNMuts=randVar.rvs(1)[0]

                    mutationsMade=0
                    while mutationsMade<actNMuts:
                        """epitope - tuple(HLA or type, seq, inds, ic50)"""
                        curEpitope=choice(epitopes)
                        """Pull the actual peptide from seq (not the insert epitope)"""
                        curEpitopeSeq=grabKmer(seq,curEpitope[2][0],k=len(curEpitope[2]))[1]
                        curEpitopeInds=grabKmerInds(seq,curEpitope[2][0],k=len(curEpitope[2]))[1]
                        if params['escapeMutations'][vaccinated]:
                            """Try all POSSIBLE single AA mutations"""
                            #print 'For PTID %d, mut %d, epitope (%s,%d-%s,%1.2f) tried variants of BT %s:' % (ptid,mutationsMade+1,curEpitope[0],curEpitope[2][0],curEpitope[1],curEpitope[3],curEpitopeSeq)
                            
                            """HAD AN ISSUE HERE WHERE curEpitopeInds WAS None?"""
                            variants=generateVariants(curEpitopeSeq,possibleAA=[possibleAA[i] for i in curEpitopeInds])
                            ic50Df=self.predictHLABinding(curEpitope[0],variants.keys())
                            """for pep,ic50 in zip(ic50Df.peptide,ic50Df.ic50):
                                print '\t%s (%1.2f)' % (pep,ic50)"""
                            btBA=ic50Df.ic50[ic50Df.peptide==curEpitopeSeq].iloc[0]
                            escapeInds=(((ic50Df.ic50-btBA) > params['escapeDelta'][vaccinated]) | (ic50Df.ic50 > 8)) & (~(ic50Df.peptide==curEpitopeSeq))
                            """Choose one of the escape variants at random"""
                            if escapeInds.sum()>=1:
                                ind=choice(ic50Df.index[escapeInds])
                                escapeVariant=ic50Df.peptide[ind]
                                posi,newAA=variants[escapeVariant]
                                sitei=curEpitopeInds[posi]
                                mutations.append((sitei,seq[sitei],newAA))
                                #print 'Chose %s to %s (%1.2f)!\n' % (curEpitopeSeq,escapeVariant,ic50Df.ic50[ind])
                                seq=mutateString(seq,sitei,newAA)
                                mutationsMade+=1
                            else:
                                """If the POSSIBLE AAs don't yield an escape then try random AAs"""
                                escapeFound=False
                                while not escapeFound:
                                    #print 'For PTID %d, mut %d(X), epitope (%s,%d-%s,%1.2f) tried variants of BT %s:' % (ptid,mutationsMade+1,curEpitope[0],curEpitope[2][0],curEpitope[1],curEpitope[3],curEpitopeSeq)
                                    variants=generateVariants(curEpitopeSeq,possibleAA=None)
                                    ic50Df=self.predictHLABinding(curEpitope[0],variants.keys())
                                    """for pep,ic50 in zip(ic50Df.peptide,ic50Df.ic50):
                                        print '\t%s (%1.2f)' % (pep,ic50)"""
                                    btBA=ic50Df.ic50[ic50Df.peptide==curEpitopeSeq].iloc[0]
                                    escapeInds=(((ic50Df.ic50-btBA) > params['escapeDelta'][vaccinated]) | (ic50Df.ic50 > 8)) & (~(ic50Df.peptide==curEpitopeSeq))
                                    if escapeInds.sum()>=1:
                                        ind=choice(ic50Df.index[escapeInds])
                                        escapeVariant=ic50Df.peptide[ind]
                                        posi,newAA=variants[escapeVariant]
                                        sitei=curEpitopeInds[posi]
                                        mutations.append((sitei,seq[sitei],newAA))
                                        print 'For PTID %d, mut %d, chose non-possible %s to %s (%1.2f) as an escape!' % (ptid,mutationsMade+1,curEpitopeSeq,escapeVariant,ic50Df.ic50[ind])
                                        seq=mutateString(seq,sitei,newAA)
                                        mutationsMade+=1
                                        escapeFound=True
                                
                        else:
                            mutationSuccess=False
                            for sitei in permutation(curEpitope[2]):
                                """Mutate site based on circulating or breakthrough strains or randomly from alphabet"""
                                #newAA=choice(alphabet)
                                tmpPossibleAA=deepcopy(possibleAA[sitei]).replace(seq[sitei],'')
                                if len(tmpPossibleAA)>0:
                                    newAA=choice(tmpPossibleAA)
                                    mutations.append((sitei,seq[sitei],newAA))
                                    seq=mutateString(seq,sitei,newAA)
                                    mutationsMade+=1
                                    mutationSuccess=True
                                    break
                            if not mutationSuccess:
                                """Optionally, could just warn and keep trying different epitopes..."""
                                raise Exception('No possible AA for epitope %s' % curEpitope)

                data['ptid'].append(str(ptid))
                data['vaccinated'].append(vaccinated)
                data['hla'].append(tmpHLA)
                data['seq'].append(seq)
                data['epitopes'].append(epitopes)
                data['mutations'].append(mutations)
                data['infected'].append(True)
                data['seqID'].append('%d_%s' % (ptid,basePtid))
                ptid+=1
            
        df=pd.DataFrame(data)
        df=df.set_index('ptid')

        """Create variables uHLA4, uHLA2 and hlaDf with 2- and 4-digit columns"""
        uHLA4=sorted(list(set(concatenate(data['hla']))))
        uHLA2=sorted(list(set([h[:-2] for h in uHLA4])))

        hlaDict={'ptid':list(df.index)}
        N=len(hlaDict['ptid'])
        for h in uHLA4+uHLA2:
            hlaDict.update({h:[False]*N})
        hlaDf=pd.DataFrame(hlaDict)
        hlaDf=hlaDf.set_index('ptid')

        for ptid in hlaDf.index:
            for h in df.hla[ptid]:
                hlaDf[h][ptid]=True
                hlaDf[h[:-2]][ptid]=True

        """Create ptidDf, seqDf, mapDf, etc."""
        self.params=params
        self.data.ptidDf=df[['vaccinated','hla','infected']]
        self.data.hlaDf=hlaDf
        self.data.seqDf=df[['seq','seqID']]
        self.data.insertSeq=sd.insertSeq
        self.data.HXB2=sd.HXB2
        self.data.mapDf=sd.mapDf
        self.data.uHLA4=uHLA4
        self.data.uHLA2=uHLA2
        self.data.N=params['n'][1]+params['n'][0]
        self.data.studyName='simulation_%s' % sd.studyName
        self.data.proteinName=sd.proteinName
        self.data.insertName='NA'
        self.data.simParams=params
        self.data.simDf=df[['mutations','epitopes']]
        self.data.date=time.asctime()
        self.data.hlaAFreq=hlaAFreq
        self.data.hlaBFreq=hlaBFreq
        self.computeDerivedData()

        """Remove HLA predictions from simulation after simulation is complete"""
        self.ba=None
        #print
    def plot(self,hxb2Range=None):
        """Plot a map of the simulated mutations and epitopes"""
      
        """Create a matrix of simulated mutations and epitopes (1/0 to indicate mutations relative to base sequence)"""
        mutMat = zeros((self.data.seqDf.shape[0],len(self.data.insertSeq),))
        epiMat = zeros((self.data.seqDf.shape[0],len(self.data.insertSeq),))
        vacInd=[]
        plaInd=[]
        for ptidi,(ptid,rec) in enumerate(self.data.simDf.iterrows()):
            if self.data.ptidDf.vaccinated[ptid]:
                vacInd.append(ptidi)
            else:
                plaInd.append(ptidi)
            for mut in rec['mutations']:
                mutMat[ptidi,mut[0]]=1
            for epi in rec['epitopes']:
                epiMat[ptidi,epi[2]]=1
        vacInd=array(vacInd)
        plaInd=array(plaInd)

        vacMut=mutMat[vacInd,:].mean(axis=0)
        plaMut=mutMat[plaInd,:].mean(axis=0)
        vacEpi=epiMat[vacInd,:].mean(axis=0)
        plaEpi=epiMat[plaInd,:].mean(axis=0)

        """Clip data to desired HXB2 range"""
        insertSeq=self.clipXVec(hxb2Range,self.data.insertSeq)
        hxb2=self.clipXVec(hxb2Range,self.data.mapDf.hxb2Pos)
        sitex=self.clipXVec(hxb2Range,self.data.mapDf.index)
        mutMat=mutMat[:,self.clipXVec(hxb2Range,returnInds=True)]
        epiMat=epiMat[:,self.clipXVec(hxb2Range,returnInds=True)]
        
        clf()
        """Plot percent mutation"""
        plot([sitex[0],sitex[-1]],[100,100],'-',color='gray')
        plot([sitex[0],sitex[-1]],[0,0],'-',color='gray')
        plot(sitex,1e2*vacMut,'-b',label='vaccine')
        plot(sitex,1e2*plaMut,'-r',label='placebo')

        """Plot percent epitope"""
        epiOffset=200
        plot([sitex[0],sitex[-1]],[-epiOffset+100,-epiOffset+100],'-',color='gray')
        plot([sitex[0],sitex[-1]],[-epiOffset,-epiOffset],'-',color='gray')
        plot(sitex,1e2*vacEpi-epiOffset,'-b')
        plot(sitex,1e2*plaEpi-epiOffset,'-r')

        mx=1e2
        yvec=arange(mutMat.shape[0])*10+mx*1.4
        
        for i,y in zip(concatenate((plaInd,vacInd)),yvec):
            ptid = self.data.seqDf.index[i]
            c = 'blue' if self.data.ptidDf.vaccinated[ptid] else 'red'
            
            for mut in self.data.simDf.mutations[i]:
                scatter(mut[0],y,c=c,alpha=0.3,s=20,marker='^',edgecolors='None')
            for ep in self.data.simDf.epitopes[i]:
                if ep[0]=='Ab':
                    for s in ep[2]:
                        plot([s-0.4,s+0.4],[y,y],'-',color='gray')
                else:
                    plot([ep[2][0],ep[2][-1]],[y,y],'-',color='gray')

        #xticks(sitex,hxb2)
        #xlim((sitex[0]-1,sitex[-1]+1))
        #xlabel('HXB2 coordinates')
        self.addHXB2Ticks(hxb2Range=hxb2Range)
        
        yticks([-epiOffset,100-epiOffset, 0,100],[0,100,0,100])
        yl=ylim()
        ylim(-epiOffset-20,yl[1]+100)

        sharedArgs=dict(xytext=(-50,0),textcoords='offset points',rotation='vertical',ha='right',va='center')
        annotate('PTID',xy=(sitex[0]-1,yvec.mean()),**sharedArgs)
        annotate('% mutated',xy=(sitex[0]-1,50),size='small',**sharedArgs)
        annotate('% epitope',xy=(sitex[0]-1,-epiOffset+50),size='small',**sharedArgs)
        title('Simulated mutations and epitopes')
        legend(loc='best')

        """
        if not saveFn is None:
            savefig(saveFn+'_sim.png',dpi=200)
        """

    def plotEpitopes(self,hxb2Range=None,sortKey='hla'):
        """
        Plot a map of the simulated epitopes
        Sort epitopes based on loc, hla, ba or freq
        """
      
        """Clip data to desired HXB2 range"""
        insertSeq=self.clipXVec(hxb2Range,self.data.insertSeq)
        hxb2=self.clipXVec(hxb2Range,self.data.mapDf.hxb2Pos)
        sitex=self.clipXVec(hxb2Range,self.data.mapDf.index)
        
        clf()
        epitopes={}
        for ptid,rec in self.data.simDf.iterrows():
            vaccinated=self.data.ptidDf.vaccinated[ptid]
            for ep in rec['epitopes']:
                #tmpEp=(ep[0],ep[1],tuple(ep[2]),ep[3])
                if not (ep[0]=='Ab' or ep[0]=='NonEpitope'):
                    if not ep in epitopes.keys():
                        epitopes.update({ep:[0,0]})
                    epitopes[ep][vaccinated]=epitopes[ep][vaccinated]+1
        
        keys=epitopes.keys()
        if sortKey=='loc':
            """Sort based on location of epitope"""
            keys=sorted(keys,key=lambda x: x[2][0])
        elif sortKey=='freq':
            """Sort based on frequency of epitope"""
            keys=sorted(keys,key=lambda x: epitopes[x][0]+epitopes[x][1])
        elif sortKey=='ba':
            """Sort based on HLA binding affinity of epitope"""
            keys=sorted(keys,key=lambda x: x[3])
        elif sortKey=='hla':
            """Sort based on HLA"""
            keys=sorted(keys,key=lambda x: x[0])

        for i,ep in enumerate(keys):
            count=epitopes[ep]
            plot(ep[2],[i]*len(ep[2]),'|-k',lw=2)
            s='%s (%1.1f)' % (ep[0].replace('_','*'),ep[3])
            annotate(s,xy=(ep[2][0],i),xytext=(-5,0),va='center',ha='right',textcoords='offset points',size='small')
            annotate('[%d,%d]' % (count[0],count[1]),xy=(ep[2][-1],i),xytext=(5,0),va='center',ha='left',textcoords='offset points',size='small')

        """Plot the vaccine sequence along x-axis"""
        ndiag=1
        for counter,aa in enumerate(self.data.insertSeq):
            annotate(s=aa,xy=(counter,-2+(ndiag-mod(counter,ndiag))),
                     ha='center',va='center',size='x-small')

        xticks(sitex,hxb2)
        yticks([])
        yl=ylim()
        ylim(yl[0]-2,yl[1]+2)
        xlim((sitex[0]-4,sitex[-1]+2))

        xlabel('HXB2 coordinates')
        sharedArgs=dict(xytext=(-50,0),textcoords='offset points',rotation='vertical',ha='right',va='center')

        title('Simulated epitopes [PLA, VAC]')

    def clipXVec(self,hxb2Range=None,vec=None,returnInds=False,clipRange=True):
        """Clip seq-axis vector based on an HXB2 coordinate range (eg [70,80])"""
        if hxb2Range is None:
            siteRange=[self.data.mapDf.index[0],self.data.mapDf.index[-1]+1]
            inds=(siteRange[0],siteRange[1]-1)
            if not vec is None:
                sitesVec=vec[siteRange[0]:siteRange[1]]
        elif clipRange:
            """Interpret the 2 element hxb2Range as a range (inclusive)"""
            hxb2Range=[str(c) for c in hxb2Range]
            siteRange=[self.data.mapDf.index[self.data.mapDf.hxb2Pos==hxb2Range[0]],self.data.mapDf.index[self.data.mapDf.hxb2Pos==hxb2Range[1]]+1]
            sitesVec=vec[siteRange[0]:siteRange[1]]
            inds=(siteRange[0],siteRange[1])
        else:
            """Interpret the hxb2Range input as specific HXB2 coordinates to slice on"""
            """HXB2 coords need to specified as the type they are in the data class mapDf"""
            #hxb2Range=[str(c) for c in hxb2Range]
            inds=[self.data.mapDf.index[self.data.mapDf.hxb2Pos==h][0] for h in hxb2Range]
            if not type(vec) is str:
                sitesVec=vec[self.data.mapDf.index[inds]]
            else:
                sitesVec=sliceString(vec,inds)
        if returnInds:
            #return arange(siteRange[0],siteRange[1])
            return inds
        else:
            return sitesVec
    def addHXB2Ticks(self,hxb2Range=None,clipRange=True):
        """Change x-axis tick marks to HXB2 coordinates"""
        
        """Labels skipping every 20 HXB2 coordinates"""
        """
        ds=20
        hxb2Labels=arange(int(self.data.mapDf.hxb2Pos[0]),int(self.data.mapDf.hxb2Pos[self.data.mapDf.shape[0]-1]),ds)-1
        hxb2Labels[0]+=1
        hxb2Labels=['%d' % s for s in hxb2Labels]
        xticks(self.data.hxb22site.site[hxb2Labels],hxb2Labels)
        """

        """Number of tick marks given by constant maxTicks"""
        maxTicks=15

        hxb2=self.clipXVec(hxb2Range,self.data.mapDf.hxb2Pos,clipRange=clipRange)
        sitex=self.clipXVec(hxb2Range,self.data.mapDf.index,clipRange=clipRange)
        
        #xt=xticks()[0]
        xt=sitex
        hxb2Labels=hxb2[xt]
        
        """I may have broke something here. I'm not sure why the HXB2 coord is required to be a str?
        I think we want them to not be null..."""
        #nani=hxb2Labels.map(lambda x: type(x) is str or type(x))
        nani=hxb2Labels.map(lambda x: type(x) is str or ~isnan(x))
        xt=xt[nani]
        xtl=hxb2Labels.dropna()
        skip=int(ceil(len(xt)/maxTicks))

        xticks(xt[::skip],xtl[::skip].tolist())
        xlim((sitex[0]-1,sitex[-1]+1))
        xlabel('HXB2 coordinate')

    def addSeq(self,hxb2Range=None,nRows=6,seq='HXB2',clipRange=True):
        """Plot the HXB2 sequence along x-axis"""
        sitex=self.clipXVec(hxb2Range,self.data.mapDf.index,clipRange=clipRange)
        consensusSeq=consensus(self.data.seqDf.seq)
        ylims=ylim()
        lineh=(ylims[1]-ylims[0])/36
        h=lineh*nRows
        counter=0
        for ix,row in self.data.hxb22site.iterrows():
            if row['site'] in sitex:
                if seq=='HXB2':
                    aa=row['hxb2aa']
                elif seq=='insert':
                    aa=self.data.insertSeq[row['site']]
                elif seq=='consensus':
                    aa=consensusSeq[row['site']]
                annotate(s=aa,xy=(row['site'],ylims[0]-h+(h/(nRows+1)*(nRows-mod(counter,nRows)-1))),
                         ha='center',va='center',size='small')
                counter+=1
        ylim((ylims[0]-lineh*(nRows+1),ylims[1]))
def generateVariants(seq,possibleAA=None):
    """Given peptide seq, return all POSSIBLE single AA mutants OR
    2 random single AA mutants per position in seq
    Return a list of variants and a dictionary of variants:(posi,newAA)"""
    variantsInfo={seq:(nan,'ORIGINAL')}
    """If we know possibleAA then return all possible variants,
    else return 2 random single-mut variant per position"""
    if not possibleAA is None:
        for posi in xrange(len(seq)):
            tmpPossibleAA=deepcopy(possibleAA[posi]).replace(seq[posi],'')
            #print ' Pos%d: %s' % (posi,tmpPossibleAA)
            for newAA in tmpPossibleAA:
                tmpVariant=mutateString(seq,posi,newAA)
                variantsInfo[tmpVariant]=(posi,newAA)
    else:
        for posi in xrange(len(seq)):
            for i in xrange(2):
                newAA=choice(AALPHABET.replace(seq[posi],''))
                tmpVariant=mutateString(seq,posi,newAA)
                variantsInfo[tmpVariant]=(posi,newAA)
    return variantsInfo

def mutateString(s,i,newChar):
    """Perform a point-mutation on string s at site i to char newChar"""
    if s[i]==newChar:
        print 'Conserved mutation %s to %s at site %d' % (s[i],newChar,i)
        raise Exception('Why')
        return s
    else:
        return s[:i] + newChar + s[i+1:]
def sliceString(s,vec):
    """Slice string s at integer inds in vec"""
    return ''.join([s[i] for i in vec])
def mutateSlice(s,vec,repl):
    """Replace s[i] with repl[j] for i in vec and j in arange(len(repl))
    Example: mutateSlice('---YY----', [3,4], 'XX') returns  '--XX----'"""
    for i,newChar in zip(vec,repl):
        s=mutateString(s,i,newChar)
    return s

def pickEpitopes(seq,btSeq,N,bindingThreshold,ba,hla,validEpitopeInds):
    """
    epitope tuple - (HLA or type, seq, inds, ic50)
    Pick N epitopes from seq using HLA-alleles in hla and IC50 in ba (dict of binding affinities)
    OR
    Pick N*9 sites at random from seq
    OR
    Pick N*9 sites based on the Brumme et al. HLA-associated sites
    """
    
    """If nEptopes is 0 then pick sites at random"""
    if N<1:
        epitopes=[('NonEpitope',seq[i],array([i]),nan) for i in permutation(validEpitopeInds)]
        inds=[e[2][0] for e in epitopes]
        return inds,epitopes

    if bindingThreshold is None:
        epitopes=[('Ab',seq[i],array([i]),nan) for i in validEpitopeInds[:(N*9)]]
        inds=[]
        for e in epitopes:
            inds+=list(e[2])
        return inds,epitopes
    else:
        """Note: grabKmerInds()[1] returns a non-gapped kmer"""
        #merInds=[tuple(grabKmerInds(seq,starti)[1]) for starti in arange(len(seq)) if not grabKmerInds(seq,starti)[1] is None]
        #mers=[grabKmer(seq,starti)[1] for starti in arange(len(seq)) if not grabKmer(seq,starti)[1] is None]
        
        merInds=[]
        mers=[]
        """Check each potential against the insert and the bt to make sure its a valid kmer
        (Prevents picking a kmer at the end of the insert that is not valid in the bt because of extra gaps at the end)"""
        for starti in arange(len(seq)):
            if not grabKmer(seq,starti)[1] is None and not grabKmer(btSeq,starti)[1] is None:
                merInds.append(tuple(grabKmerInds(seq,starti)[1]))
                mers.append(grabKmer(seq,starti)[1])
        
        potentialEpitopes=[]
        for h in hla:
            for pep,ind in zip(mers,merInds):
                """Require that at least 1 residue can tolerate a mutation"""
                if any([(meri in validEpitopeInds) for meri in ind]):
                    tmpba=ba[(h,pep)]
                    if tmpba<bindingThreshold:
                        potentialEpitopes.append((h,pep,ind,tmpba))

        if len(potentialEpitopes)>=N:
            #return unique epitopes...
            epitopes=[]
            inds=[]
            for i in arange(N):
                tmp=choice(potentialEpitopes)
                inds+=list(tmp[2])
                epitopes.append(tmp)
                potentialEpitopes.remove(tmp)
            return inds,epitopes
            #this line may return non-unique epitopes        
            #return [choice(potentialEpitopes)[2] for i in range(N)]
        else:
            return [],[]
    

"""
descriptionStr='mutations_per_epitope=%d|N_epitopes=%d|vacThreshold=%d|mutation_prob=%1.2f' % (mutations_per_epitope,N_epitopes,vac_epitope_threshold,mutation_prob)
"""

def simSDFromLANL(dataPath,protein,year,hlaFreq,clade=None,country=None):
    """Starts with an empty sieveData object and fills it in preparation for a simulation requiring these fields only:
            insertSeq
            hlaDf (or hlaFreq specified)
            ptidDf - will limit to placebo sequences (column: ~vaccinated)
            HXB2
            mapDf
            studyName
            proteinName
        Try: simSD=simSDFromLANL(protein='gag',year=2012,clade='C',hlaFreq=loadSouthAfricaHLA())
    """
    sd=sieveData()

    sd.proteinName=protein
    sd.studyName='LANL_Clade%s_%d' % (clade,year)

    lanlDf=fasta2df(dataPath + 'HIV_alignments/HIV1_FLT_%s_%s_PRO.fasta' %(year,protein))

    sd.HXB2=lanlDf.seq[lanlDf.name=='HXB2_LAI_IIIB_BRU'].tolist()[0]
    if not clade is None:
        lanlDf=lanlDf.ix[lanlDf.clade==clade]
    if not country is None:
        lanlDf=lanlDf.ix[lanlDf.country==country]
    for ind in lanlDf.index:
        lanlDf.seq[ind]=re.sub('[%s]' % BADAA,'-',lanlDf.seq[ind])
    lanlDf.seq=padAlignment(lanlDf.seq)

    """Remove any sequences that have a gap at the positions that only have 1, 2 or 3 gaps.
    This increases the number of sites that can be analyzed/simulated without reducing the number of sequences by too much
    (With <=3 it should still be 396 clade C ZA in 2012, while for <=100 its 277 valid seqs)"""
    """smat=align2mat(lanlDf.seq)
    gapSiteInd=(sum(smat=='-',axis=0)>0) & (sum(smat=='-',axis=0)<=100)
    gapSeqInd=any(smat[:,gapSiteInd]=='-',axis=1)
    lanlDf=lanlDf.ix[lanlDf.index[~gapSeqInd]]"""

    lanlDf['vaccinated']=zeros(lanlDf.shape[0],dtype=bool)
    lanlDf['infected']=ones(lanlDf.shape[0],dtype=bool)
    lanlDf['ptid']=['%s' % s for s in xrange(lanlDf.shape[0])]

    lanlDf=lanlDf.set_index('ptid')[['seq','vaccinated','infected']]
    sd.seqDf=lanlDf[['seq']]
    sd.ptidDf=lanlDf[['vaccinated','infected']]
    
    """Use consensus or mindist for teh insert sequence"""
    #sd.insertSeq=consensus(sd.seqDf.seq)
    sd.insertSeq=identifyMindist(sd.seqDf.seq,ignoreGaps=False)
    
    """Set up mapDf"""
    sd.mapDf=pd.DataFrame(zeros((len(sd.insertSeq),1),dtype=int32),columns=['posNum'])
    sd.mapDf.posNum=arange(len(sd.insertSeq))
    sd.mapDf['hxb2AA']=array([aa for aa in sd.HXB2],dtype='S1')
    sd.mapDf['hxb2Pos']=array(['' for aa in sd.HXB2],dtype=object)
    sd.mapDf=sd.mapDf.set_index('posNum')
    tmpCount=0
    lettLookup={n:a for n,a in zip(arange(26),'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.lower())}
    lettCount=0
    for sitei in sd.mapDf.index:
        if sd.mapDf.hxb2AA[sitei]=='-':
            sd.mapDf.hxb2Pos[sitei]='%d%s' % (tmpCount,lettLookup[lettCount])
            lettCount+=1
        else:
            tmpCount+=1
            lettCount=0
            sd.mapDf.hxb2Pos[sitei]='%d' % (tmpCount)

    sd.hlaFreq=hlaFreq
    uHLA=sd.hlaFreq['A'].keys()+sd.hlaFreq['B'].keys()
    sd.uHLA4=[h[:6] for h in uHLA]
    sd.uHLA2=[h[:4] for h in uHLA]
    sd.lookupDf=None
    sd.hxb22site=None
    return sd