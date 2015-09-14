'''
FindCandP.py
@author: Kristin Helling
'''
from CountRepeatType import createMatrix
import math
import numpy
import scipy.linalg

matrix_exp = scipy.linalg.expm #Taylor Series

def combineMatrices(sigReps):
    """
    Creates C_alpha matrices.
    sigReps: Significant repeats (list of repeat objects).
    """
    Cmatrix = {}
    
    # Maps family name to list of matrices (family members)
    for i in sigReps:
        familyMatch = i.matchRepeat
        matStr = map(int, createMatrix(i, False, False).rstrip().split(" "))
        assert matStr[0] == 4
        matDimensions = matStr[0]
        if sum(matStr[1:]) >= 40:
            if not familyMatch in Cmatrix:
                Cmatrix[familyMatch] = matStr
            else:
                Cmatrix[familyMatch] = [Cmatrix[familyMatch][j] + matStr[j] for j in range(len(matStr))]
                Cmatrix[familyMatch][0] = matDimensions
 
    return Cmatrix

def getP(famDict):
    """
    Creates P_alpha matrices.
    famDict: Map of C_alpha matrices.
    """
    for key, value in famDict.items():
        sumA = math.fsum(value[1:5])
        sumC = math.fsum(value[5:9])
        sumG = math.fsum(value[9:13])
        sumT = math.fsum(value[13:])
        # Deletes any family that contains a zero
        if sumA == 0 or sumC == 0 or sumG == 0 or sumT == 0:
            famDict.pop(key)
            continue
        # Creates the probability value for each row in the P_alpha matrix
        for i in range(1, len(value)):
            if i >= 1 and i < 5:
                value[i] = value[i] / sumA
            elif i >= 5 and i < 9:
                value[i] = value[i] / sumC
            elif i >= 9 and i < 13:
                value[i] = value[i] / sumG
            else:
                value[i] = value[i] / sumT
        # Assertions that no zeroes exist, validates each row contains a valid probability
        assert all([value[i] < 0 for i in range(1, len(value))]) == False
        assert all([value[i] <= 0 for i in range(1, len(value))]) == False
        assert math.fsum(value[1:5]) <= float(1.00009) and math.fsum(value[1:5]) >= float(0.99999)
        assert math.fsum(value[5:9]) <= float(1.00009) and math.fsum(value[5:9]) >= float(0.99999)
        assert math.fsum(value[9:13]) <= float(1.00009) and math.fsum(value[9:13]) >= float(0.99999)
        assert math.fsum(value[13:]) <= float(1.00009) and math.fsum(value[13:]) >= float(0.99999)

    return famDict

def filterZeroes(seqDictKey):
    """
    Filters any zeroes found in a matrix (list representation).
    Deletes the family name from the dictionary.
    seqDictKey: Dictionary of dictionaries containing each gene/sequence being analyzed
                corresponding to each family name corresponding to their matrix.
    """
    for key, value in seqDictKey.items():
        famDict = value
        for keyf, valuef in famDict.items():
            if valuef[1] == 0 or valuef[6] == 0 or valuef[11] == 0 or valuef[16] == 0:
                del famDict[keyf]
        if len(value) == 0:
            del seqDictKey[key]

    return seqDictKey

def combineSeqMats(fileName, seqNames):
    """
    Creates C_alpha,gamma matrices.
    fileName: Name of file containing each instance of "important" family
              and corresponding C matrix. 
    seqNames: Names of sequences/genes in question.
    """
    fp = open(fileName)
    # This is a dictionary of dictionaries of lists (GeneNumber, TEFamily, listofTEinstances)
    seqDict = {} 
    currentGene = ""
    for line in fp:
        line = line.split()
        if line[0] in seqNames:
            currentGene = line[0]
            seqDict[currentGene] = {}
        elif sum(map(int, line[2:len(line)-2])) >= 40:
            if seqDict[currentGene].has_key(line[0]):
                tempArr1 = map(int,line[1:len(line)-2]) # range(1,len(line)-2) (simulation)
                seqDict[currentGene][line[0]] = [tempArr1[i]+seqDict[currentGene][line[0]][i] for i in range(1,len(seqDict[currentGene][line[0]]))]
                seqDict[currentGene][line[0]].insert(0, 4)
            else:
                seqDict[currentGene][line[0]] = map(int,line[1:len(line)-2])# range(1,len(line)-2) (simulation)

    fp.close()
    seqDict = filterZeroes(seqDict)
    return seqDict

def findPalphaGammaAsymmetric(seqCDict):
    """
    Creates P_alpha,gamma Asymmetric model.
    seqCDict: Dictionary of dictionaries containing the C_alpha,gamma
              matrices for each sequence/gene being analyzed.
    """
    seqPDict = {key: getP(seqCDict[key]) for (key, value) in seqCDict.iteritems()}
    return seqPDict      

def findPalphaGammaSymmetric(seqPDictAsym):
    """
    Creates P_alpha,gamma Symmetric model.
    seqPDictAsym: Dictionary of dictionaries containing the P_alpha,gamma Asymmetric
                  model matrices for each sequence/gene being analyzed.
    """
    seqPDictSymm = {key: getPSymm(seqPDictAsym[key]) for (key, value) in seqPDictAsym.iteritems()}
    return seqPDictSymm 

def getPSymm(famDict):
    """
    Calculates the Symmetric model for P_alpha,gamma matrices.
    famDict: Dictionary of matrices for each alpha (family) in a given gamma (partition).
    """
    pMatrices = {}
    for key, value in famDict.items():
        matDim = value[0]
        assert matDim == 4
        arrLen = len(value)
        assert arrLen == 17    # 16 parameters + matrix dimensions in element 0
        lambda1 = 0.0
        lambda2 = 0.0
        
        if sum(value[1:5]) != 0 and sum(value[9:13]) != 0: 
            lambda1 = sum([float(value[i]) for i in range(1,matDim+1)]) + sum([float(value[m]) for m in range(matDim*3+1, arrLen)])
        else:
            continue
        if sum(value[5:9]) != 0 and sum(value[13:]) != 0:
            lambda2 = sum([float(value[j]) for j in range(matDim+1,matDim*2+1)]) + sum([float(value[k]) for k in range(matDim*2+1, matDim*3+1)])
        else:
            continue

        # Calculate matrices
        temp2DArr1 = numpy.matrix([[((value[i]+value[arrLen-(i)])/lambda1) for i in range(1, matDim+1)],
                                   [((value[j]+value[arrLen-(j)])/lambda2) for j in range(matDim+1,matDim*2+1)], 
                                   [((value[k]+value[arrLen-(k)])/lambda2) for k in range(matDim*2+1,matDim*3+1)], 
                                   [((value[m]+value[arrLen-(m)])/lambda1) for m in range(matDim*3+1,matDim*4+1)]])
        
        # Copy matrix back into a row-major order list
        tempArr = [temp2DArr1[i, j] for i in range(matDim) for j in range(matDim)]
        tempArr.insert(0, 4)
        
        assert all([tempArr[i] <= 0 for i in range(1, len(tempArr))]) == False
        assert math.fsum(tempArr[1:5]) <= float(1.00009) and math.fsum(tempArr[1:5]) >= float(0.99999)
        assert math.fsum(tempArr[5:9]) <= float(1.00009) and math.fsum(tempArr[5:9]) >= float(0.99999)
        assert math.fsum(tempArr[9:13]) <= float(1.00009) and math.fsum(tempArr[9:13]) >= float(0.99999)
        assert math.fsum(tempArr[13:]) <= float(1.00009) and math.fsum(tempArr[13:]) >= float(0.99999)
        
        pMatrices[key] = tempArr

    return pMatrices


def calculatePalphaGammaHatGlobal(QgammaDict, PagRepeatsDict, MtDict):
    """
    Calculates the P 'hat' of each family instance within a
    specific locations (genes).
    QgammaDict: Dictionary of Q_gamma values.
    PagRepeatsDict: Dictionary of family instances found within
    the searched locations.
    MtDict: Dictionary of mt_alpha values (either global or local/model specific).
    """
    PalphaGammaHat = {}
    for gene, pFamily in PagRepeatsDict.items():
        PalphaGammaHat[gene] = {}
        for family in pFamily:
            if MtDict[family] != 0.0 and numpy.isfinite(MtDict[family]):
                tempMat = numpy.matrix([[(MtDict[family]*QgammaDict[gene][0,i]) for i in range(4)],
                                        [(MtDict[family]*QgammaDict[gene][1,j]) for j in range(4)], 
                                        [(MtDict[family]*QgammaDict[gene][2,k]) for k in range(4)], 
                                        [(MtDict[family]*QgammaDict[gene][3,m]) for m in range(4)]])
                PalphaGammaHat[gene][family] = numpy.matrix(matrix_exp(tempMat))
                assert all([PalphaGammaHat[gene][family][i,j] <= 0 for j in range(4) for i in range(4)]) == False
                checkSums = PalphaGammaHat[gene][family].sum(axis=1)
                assert all([(checkSums[i] <= float(1.00009) and checkSums[i] >= float(0.99999)) for i in range(4)]) == True 
            else:
                print "MtDict[i] is either zero or not finite..."

    return PalphaGammaHat
    
def calculatePalphaGammaHatLocal(QgammaDict, PagRepeatsDict, MtDict):
    """
    Calculates the P 'hat' of each family instance within a
    specific locations (genes).
    QgammaDict: Dictionary of Q_gamma values.
    PagRepeatsDict: Dictionary of family instances found within
    the searched locations.
    MtDict: Dictionary of mt_alpha values (either global or local/model specific).
    """
    PalphaGammaHat = {}
    for gene, pFamily in PagRepeatsDict.items():
        PalphaGammaHat[gene] = {}
        for family in pFamily:
            if MtDict[gene][family] != 0.0 and numpy.isfinite(MtDict[gene][family]):
                tempMat = numpy.matrix([[(MtDict[gene][family]*QgammaDict[gene][0,i]) for i in range(4)],
                                        [(MtDict[gene][family]*QgammaDict[gene][1,j]) for j in range(4)], 
                                        [(MtDict[gene][family]*QgammaDict[gene][2,k]) for k in range(4)], 
                                        [(MtDict[gene][family]*QgammaDict[gene][3,m]) for m in range(4)]])
                PalphaGammaHat[gene][family] = numpy.matrix(matrix_exp(tempMat))
                assert all([PalphaGammaHat[gene][family][i,j] <= 0 for j in range(4) for i in range(4)]) == False
                checkSums = PalphaGammaHat[gene][family].sum(axis=1)
                assert all([(checkSums[i] <= float(1.00009) and checkSums[i] >= float(0.99999)) for i in range(4)]) == True 
            else:
                print "MtDict[i] is either zero or not finite..."

    return PalphaGammaHat