'''
printAndParseFiles.py
@author: Kristin Helling
'''
import re
import numpy

def printListToFile(listDictionary, fileName, fp = ""):
    """
    Prints a dictionary of lists to a file.
    listDictionary: Dictionary of lists.
    fileName: Name of output file.
    fp: File pointer, if one already exists.
    """
    needClose = False
    if fp == "":
        fp = open(fileName, 'w')
        needClose = True
    for key, value in listDictionary.items():
        line = key + " " + " ".join([str(i) for i in value])
        fp.write(line + "\n")
    if needClose == True:
        fp.close()
    
def printDictionaryListToFile(dictionaryList, fileName):
    """
    Prints a dictionary of dictionaries of lists to a file.
    dictionaryList: Dictionary of dictionaries of lists.
    fileName: Name of output file.
    """
    fp = open(fileName, 'w')
    for key, value in dictionaryList.items():
        fp.write(key + "\n")
        printListToFile(value, fileName, fp)
    fp.close()
    
def printMatrixToFile(matrixDictionary, fileName):
    """
    Prints a dictionary of matrices to a file.
    matrixDictionary: Dictionary of matrices.
    fileName: Name of output file.
    """
    fp = open(fileName, 'w')
    for key, value in matrixDictionary.items():
        line = key + " " + " ".join([str(value[row,col]) for row in range(len(value)) for col in range(len(value))])
        fp.write(line + "\n")
    fp.close()

def printDictionaryMatrixToFile(dictionaryMatrix, fileName):
    """
    Prints a dictionary of dictionaries of matrices to a file.
    dictionaryMatrix: Dictionary of dictionaries of matrices.
    fileName: Name of output file.
    """
    fp = open(fileName, 'w')
    for sequence, family in dictionaryMatrix.items():
        fp.write(sequence + "\n")
        for key, value in family.items():
            line = key + " " + " ".join([str(value[row,col]) for row in range(len(value)) for col in range(len(value))])
            fp.write(line + "\n")
    fp.close()

def printValuesToFile(valueDictionary, fileName, fp = ""):
    """
    Prints a dictionary of values to a file.
    valueDictionary: Dictionary of values (scalar).
    fileName: Name of output file.
    """
    needClose = False
    if fp == "":
        fp = open(fileName, 'w')
        needClose = True
    for key, value in valueDictionary.items():
        line = key + " " + str(value) + "\n"
        fp.write(line)
    if needClose == True:
        fp.close()

def printDictionaryofValuesToFile(valueDictionary, fileName):
    """
    Prints a dictionary of dictionaries of values to a file.
    valueDictionary: Dictionary of values (scalar).
    fileName: Name of output file.
    """
    fp = open(fileName, 'w')
    for key, value in valueDictionary.items():
        fp.write(key + "\n")
        printValuesToFile(value, fileName, fp)
    fp.close()
    
def compareAndPrintBICs(BICA, BICS, fileName):
    """
    Compares and outputs the determination of 
    Asymmetric or Symmetric to the given gamma.
    If the BICS contains other keys the BICA does
    not contain, return the key as 'Symmetric'.
    BICA: Dictionary of values with the Asymmetric model calculation.
    BICS: Dictionary of values with the Symmetric model calculation.
    fileName: Output file name.
    """
    fp = open(fileName, 'w')
    genes = len(BICA)
    count = 0
    for key, value in BICA.items():
        if BICS.has_key(key):
            if BICS[key] > BICA[key]:
                fp.write(key + " Asymmetric\n")
                count = count + 1
            elif BICS[key] < BICA[key]:
                fp.write(key + " Symmetric\n")
        else:
            fp.write(key + " Asymmetric\n")
            count = count + 1
    for key, value in BICS.items():
        if not BICA.has_key(key):
            fp.write(key + " Symmetric\n")
            genes = genes + 1
            
    fp.close()
    print "COUNT: " + str(count) + " are Asymmetric out of " + str(genes)

def parsePfamilyFile(fileName):
    """
    Parses P_alpha file in order to save time after
    first calculation of C_alpha and P_alpha.
    fileName: File containing P_alpha data.
    """
    fp = open(fileName)
    PfamilyMats = {line.split()[0]: map(float,line.split()[1:]) for line in fp}
    fp.close()
    return PfamilyMats

def SigTEs(allRepObjs, repeatFamilyFile):
    """
    Creates list of significant Repeat objects.
    allRepObjs: Complete list of Repeat objects.
    repeatFamilyFile: File containing 'important' families.
    """
    # opens and reads repeat list (contains family names)
    with open(repeatFamilyFile) as repClassNames:
        reptList = [line.strip() for line in repClassNames]

    # keeps the alignments found in repeat family file
    keepAlign = [i for i in allRepObjs if i.matchRepeat in reptList]  

    return keepAlign

def findSeqNames(fileName):
    """
    Reads in a file with all of the sequence/gene names.
    Returns a list of names.
    fileName: Name of file containing sequence/gene names.
    """
    fp = open(fileName)
    seqNames = [line.split()[0] for line in fp]
    fp.close()
    return seqNames

def parsePalphaGamma(fileName, beginSeqName):
    """
    Parses the P_alpha,gamma file for repeats within a sequence.
    Returns a dictionary with sequence query names and corresponding
    repeats.
    fileName: file containing P_alpha,gamma data.
    beginSeqName: beginning of each sequence name for 
                  the given query data (e.g. "NM_").
    """
    fp = open(fileName)
    repeatDict = {}
    for line in fp:
        if re.search("^" + beginSeqName, line):
            sequence = line.strip()
            repeatDict[sequence] = []
        else:
            repeatDict[sequence].append(line.split()[0])
    fp.close()

    return repeatDict

def parseQgammaFile(fileName):
    """
    Parses the Q_gamma values found within the
    file given.
    fileName: Q_gamma file.
    """
    fp = open(fileName)
    Qgamma = {}
    for line in fp:
        line = line.split()
        Qgamma[line[0]] = numpy.matrix([[float(line[i]) for i in range(1,5)],
                                        [float(line[j]) for j in range(5,9)], 
                                        [float(line[k]) for k in range(9,13)], 
                                        [float(line[m]) for m in range(13,17)]])
    fp.close()
    return Qgamma  

def parseCgammaFile(fileName):
    """
    Parses the C_gamma values found within the
    file given.
    fileName: C_gamma file.
    """
    fp = open(fileName)
    Cgamma = {}
    for line in fp:
        line = line.split()
        Cgamma[line[0]] = numpy.matrix([[int(line[i]) for i in range(1,5)],
                                        [int(line[j]) for j in range(5,9)], 
                                        [int(line[k]) for k in range(9,13)], 
                                        [int(line[m]) for m in range(13,17)]])
    fp.close()
    return Cgamma  

    
def parseValueFile(fileName):
    """
    Parses the values found within the
    file given.
    fileName: Either mt_alpha or L_gamma file.
    """
    fp = open(fileName)
    valDict = {line.split()[0]: float(line.split()[1]) for line in fp}
    fp.close()
    return valDict

def parseCalphaGammaFile(fileName, beginSeqName):
    """
    Parses the C_alpha,gamma values found within the
    file given.
    fileName: C_alpha,gamma file.
    beginSeqName: beginning of each sequence name for 
                  the given query data (e.g. "NM_").
    """
    fp = open(fileName)
    CalphaGamma = {}
    for line in fp:
        if re.search("^" + beginSeqName, line):
            repeatDict = {}
            CalphaGamma[line.strip()] = repeatDict
        else:
            line = line.split()
            repeatDict[line[0]] = numpy.matrix([[int(line[i]) for i in range(2,6)],
                                                [int(line[j]) for j in range(6,10)], 
                                                [int(line[k]) for k in range(10,14)], 
                                                [int(line[m]) for m in range(14,18)]])
    fp.close()
    return CalphaGamma  

def parseMtalphaGammaFile(fileName, beginSeqName):
    fp = open(fileName)
    MtalphaGamma = {}
    for line in fp:
        if re.search("^" + beginSeqName, line):
            repeatDict = {}
            MtalphaGamma[line.strip()] = repeatDict
        else:
            line = line.split()
            repeatDict[line[0]] = float(line[1].strip())
    fp.close()
    return MtalphaGamma

def parseBICfile(fileName):
    fp = open(fileName)
    valDict = {line.split()[0]: line.split()[1] for line in fp}
    fp.close()
    return valDict