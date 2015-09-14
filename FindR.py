'''
FindR.py
@author: Kristin Helling
'''
import numpy
import scipy.linalg
eps = numpy.finfo(float).eps
import scipy.linalg._matfuncs_inv_ssq

def calculateR(Pmatrix, chr = "", gene = "", symmetry = "", directory = ""):
    """
    Creates dictionary of lists representing R_alpha matrices.
    Pmatrix: Dictionary of lists of P_alpha matrices.
    chr: Chromosome number represented as a string (used for R_alpha,gamma calculation).
    gene: Name of the current gene (used for R_alpha,gamma calculation).
    symmetry: String representation of which model currently being calculated 
              ("Asymmetric" or "Symmetric").
    directory: Name of the directory for filtered family output files for this
               calculation. 
    """
    RmatDict = {}
    fpDiag = open(directory + "DiagonalAssertionFailureKeys" + chr + symmetry + ".txt", "a")
    fpOffDiag = open(directory + "OffDiagonalAssertionFailureKeys" + chr + symmetry + ".txt", "a")
    if not gene == "":
        fpDiag.write(gene + "\n")
        fpOffDiag.write(gene + "\n")

    for key, value in Pmatrix.items():
        # Creates the matrix structure based on the area dimensions
        matrixArea = int(value[0])
        assert matrixArea == 4

        Pmat = numpy.matrix([[value[i] for i in range(1,1+matrixArea)],
                             [value[j] for j in range(1+matrixArea,1+matrixArea*2)], 
                             [value[k] for k in range(1+matrixArea*2,1+matrixArea*3)], 
                             [value[m] for m in range(1+matrixArea*3,1+matrixArea*4)]])
        
        # similar to scipy, used to identify problem matrices
        F = scipy.linalg.logm(Pmat)
        errtol = 1000*eps
        errest = scipy.linalg.norm(scipy.linalg.expm(F)-Pmat,1) / scipy.linalg.norm(Pmat,1)
        if numpy.linalg.det(Pmat) == 0:
            print "Determinant is zero for Chromosome " + chr + ", " + symmetry + ", " + gene + ", " + key
        w, v = numpy.linalg.eig(Pmat)
        if 0.0 in w:
            print "Zero in eigenvalues for Chromosome " + chr + ", " + symmetry + ", " + gene + ", " + key 
        if not numpy.isfinite(errest):
            print "Is not finite for Chromosome " + chr + ", " + symmetry + ", " + gene + ", " + key 
        if errest >= errtol:
            print "logm result may be inaccurate, approximate error =" + str(errest) + " for Chromosome " + chr + ", " + symmetry + ", " + gene + ", " + key

        Rmatrix = scipy.linalg.logm(Pmat)

        Rmatrix = numpy.matrix([ [x.real for x in row] for row in Rmatrix ])
        # Assertions that all calculations return valid numbers
        assert all([isinstance(Rmatrix[i, j], complex) for j in range(4) for i in range(4)]) == False	
        if all([Rmatrix[i, i] < 0 for i in range(4)]) == False:
            fpDiag.write(key + "\n")
        elif all([Rmatrix[i, j] > (-1.0e-5) for j in range(4) for i in range(4) if i != j]) == False:
            fpOffDiag.write(key + "\n")
        else:
            RmatDict[key] = Rmatrix

    fpDiag.close()
    fpOffDiag.close() 

    return RmatDict

def calculateRag(PagDictionary, chr, symmetry, directory):
    """
    Creates dictionary of dictionaries for R_alpha,gamma matrices.
    Pagmatrix: Dictionary of dictionaries representing P_alpha,gamma matrices.
    chr: Chromosome number represented as a string.
    symmetry: String representation of which model currently being calculated 
              ("Asymmetric" or "Symmetric").
    directory: Name of the directory for filtered family output files for this
               calculation.
    """
    RagDictionary = {key: calculateR(value, chr, key, symmetry, directory) for (key, value) in PagDictionary.items()}
    return RagDictionary
