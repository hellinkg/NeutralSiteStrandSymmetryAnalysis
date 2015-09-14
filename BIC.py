'''
BIC.py
@author: Kristin Helling
'''
import numpy
    
def findBIC(PagDict, LDict, AorS):
    """
    Calculates the BIC value for each gamma,
    depending upon the model being used.
    PagDict: P_alpha,gamma hat dictionary.
    LDict: L_gamma dictionary.
    AorS: Character indicating Asymmetric ("A") or
          Symmetric ("S") model.
    """
    BIC = {}
    if AorS == "A":
        f = 11.0
        BIC = {gene: (-2.0*value) + ((f+len(PagDict[gene])) * numpy.log((16.0*len(PagDict[gene])))) for gene,value in LDict.items()}
    elif AorS == "S":
        f = 5.0
        BIC = {gene: (-2.0*value) + ((f+len(PagDict[gene])) * numpy.log((16.0*len(PagDict[gene])))) for gene,value in LDict.items()}
    return BIC