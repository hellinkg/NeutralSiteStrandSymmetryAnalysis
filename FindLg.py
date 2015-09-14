'''
FindLg.py
@author: Kristin Helling
'''
import numpy
import math

def findLgamma(CalphaGammaDict, PalphaGammaHatDict):
    """
    Calculates L_gamma using the previously calculated
    C_alpha,gamma and P_alpha,gamma 'hat' dictionaries.
    CgammaDict: Dictionary of dictionaries representing C_alpha,gamma values.
    PgammaDict: Dictionary of dictionaries representing P_alpha,gamma 'hat' values.
    """
    LgammaDict = {}
    for gene, value in PalphaGammaHatDict.items():
        LgammaValue = 0.0
        for family in value:
            if family in CalphaGammaDict[gene]:
                if all(PalphaGammaHatDict[gene][family][row,col] != 0.0 for row in range(4) for col in range(4)):
                    LgammaValue += numpy.sum(numpy.matrix([[CalphaGammaDict[gene][family][0,i]*numpy.log(PalphaGammaHatDict[gene][family][0,i]) for i in range(4)],
                                                           [CalphaGammaDict[gene][family][1,j]*numpy.log(PalphaGammaHatDict[gene][family][1,j]) for j in range(4)], 
                                                           [CalphaGammaDict[gene][family][2,k]*numpy.log(PalphaGammaHatDict[gene][family][2,k]) for k in range(4)], 
                                                           [CalphaGammaDict[gene][family][3,m]*numpy.log(PalphaGammaHatDict[gene][family][3,m]) for m in range(4)]]))
        if LgammaValue != 0.0 and not math.isnan(LgammaValue):
            LgammaDict[gene] = LgammaValue
    return LgammaDict