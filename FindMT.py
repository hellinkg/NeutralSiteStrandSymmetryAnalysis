'''
FindMT.py
@author: Kristin Helling
'''
import numpy

def calculateMt(Rmatrix, matSize):
    """
    Calculates the mt_alpha values.
    Rmatrix: R_alpha matrices.
    matSize: Dimensions of the matrices.
    """
    assert matSize == 4
    # Uses numpy to sum the diagonal of the R matrices
    # in order to calculate mt_alpha
    matDict = {key: (numpy.trace(value) / matSize) * -1 for (key, value) in Rmatrix.items()}
    return matDict    

def calculateMtalphaGamma(RagDictionary, matSize):
    """
    Calculates the mt_alpha,gamma values.
    Rmatrix: Dictionary of dictionaries representing R_alpha,gamma matrices.
    matSize: Dimensions of the matrices.
    """
    MtagDictionary = {key: calculateMt(value, matSize) for (key, value) in RagDictionary.items()}
    return MtagDictionary