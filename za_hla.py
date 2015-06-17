"""
Module for parsing HLA frequencies from a patient data file downloaded from
one of the LANL epitope mapping data sets.

http://www.hiv.lanl.gov/content/immunology/hlatem/study4/index.html

"""

import pandas as pd
from objhist import objhist

__all__ = ['loadSouthAfricanHLA']

def loadSouthAfricaHLA(filePath):
    """Returns a dict of HLA frequencies that can be used
    to simulate patient HLAs in pysieve/simulation.py"""
    #DATA_PATH + 'HLATEM/study4/patients.csv'
    rawDf = pd.read_csv(filePath).dropna()
    return {'A':_processHLAs(rawDf,['HLA_A1','HLA_A2']),
            'B':_processHLAs(rawDf,['HLA_B1','HLA_B2'])}
    
def _processHLAs(rawDf, columns, hlaSubs = None):
    """This parsing works for 'HLATEM/study4/patients.csv'"""
    if hlaSubs is None:
        hlaSubs=[('A*11','A*1101'),
                 ('A*28','A*6801'),
                 ('A*26','A*2601'),
                 ('A*202','A*0202'),
                 ('A*39',None),
                 ('A*0702',None),
                 ('B*27','B*2702'),
                 ('B*35','B*3501'),
                 ('B*54','B*5401'),
                 ('B*71',None),
                 ('B*72',None),
                 ('B*0202',None)]
    def replaceKey(a,oldk,newk):
        if a.has_key(oldk):
            if not newk is None:
                tmp = a[oldk]
                try:
                    a[newk] += tmp
                except KeyError:
                    a[newk] = tmp
            dropped = a.pop(oldk)

    a = objhist(rawDf[columns[0]])
    for c in columns[1:]:
        a.add(rawDf[c])

    for k in a.keys():
        if not isinstance(k,basestring):
            dropped = a.pop(k)
            continue
        
        if k.find('/')>=0:
            newk = k.split('/')[0]
            replaceKey(a,k,newk)
            k = newk
        if k.find('\\')>=0:
            newk=k.split('\\')[0]
            replaceKey(a,k,newk)
            k=newk
        if '?' in k:
            newk = k.replace('?','')
            replaceKey(a,k,newk)
            k = newk

    for oldk,newk in hlaSubs:
        replaceKey(a,oldk,newk)

    """Go through HLAs that are < 4-digit and replace with most common match"""
    replacements = []
    for k in a.keys():
        if len(k) < 6:
            mx = 0
            mxk = None
            for newk in a.keys():
                if len(newk) == 6 and newk[:len(k)] == k:
                    mx = a[newk]
                    mxk = newk
            if mxk is None:
                print 'No match for %s!' % k
            else:
                replacements.append((k,mxk))
        if len(k) > 6:
            replacements.append((k,k[:6]))
    for oldk,newk in replacements:
        replaceKey(a,oldk,newk)
    return a.freq()