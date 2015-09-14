'''
FindQ.py
@author: Kristin Helling
'''
from collections import defaultdict
import numpy
import math

def calculateQvalues(mtDictionary, Rmatrix, matSize):
    """
    Calculates Q_alpha values.
    mtDictionary: Dictionary of mt_alpha values.
    Rmatrix: Dictionary of R_alpha values.
    matSize: Dimensions of R_alpha values.
    """
    assert matSize == 4
    Qmatrix = {}
    
    # Divides R_alpha values by the corresponding mt_alpha value
    for key, value in Rmatrix.items():
        if key in mtDictionary and numpy.isfinite(mtDictionary[key]):
            Qmatrix[key] = numpy.divide(value, mtDictionary[key])

        # Assertion that ensures each row sums to approximately zero
        assertMat = Qmatrix[key].sum(axis=1)
        assert all([(assertMat[i] <= float(.00009) and assertMat[i] >= float(-0.00009)) for i in range(matSize)])
    return Qmatrix

def calculateQag(RagDictionary, mtDictionary, matSize):
    """
    Creates Q_alpha,gamma matrices.
    RagDictionary: Dictionary of dictionaries representing R_alpha,gamma matrices.
    mtDictionary: Dictionary of mt_alpha,gamma values.
    """
    QagDictionary = {key: calculateQvalues(mtDictionary[key], value, matSize) for (key,value) in RagDictionary.items()}
                
    return QagDictionary
    
def findQgamma(QagDictionary):
    """
    Creates Q_gamma matrices.
    QagDictionary: Dictionary of dictionaries representing Q_alpha,gamma matrices.
    """
    QgammaDict = defaultdict(float)
    for sequence, repeats in QagDictionary.items():
        if len(repeats) != 0:
            for key, value in repeats.items():
                QgammaDict[sequence] += value 
            QgammaDict[sequence] = numpy.divide(QgammaDict[sequence],len(repeats.items()))

    return QgammaDict