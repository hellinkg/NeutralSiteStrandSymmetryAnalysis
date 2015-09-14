'''
RunAllCalcs.py
@author: Kristin Helling
'''
# Helper Library/File
import sys
import itertools
import copy
from RMfileReader import *

# All Calculations
import FindCandP
import FindR
import FindMT
import FindQ
import FindLg

# BIC Test
import BIC

# File I/O
import printAndParseFiles

def RunGlobalAlphaTest(beginSeqNum, endSeqNum, firstTime, significantRepeatList, organism, prmPrefix, prmDirectory = 'prm_files/'):
    """
    Runs all Global Alpha test calculations.
    beginSeqNum: Beginning chromosome range.
    endSeqNum: Ending chromosome range.
    firstTime: Boolean if first time for all calculations.
    significantRepeatList: List of significant transposable element family instances.
    organism: Organism's name.
    prmDirectory: Directory where all of the .prm files are located.
    """
    print ("GLOBAL ALPHA")
    # to find global alpha values for C and P
    if firstTime == True:
        AllRepObjs = list(itertools.chain.from_iterable(([unpickleRepeats(prmPrefix + "%d" % chrNum, prmDirectory) for chrNum in range(beginSeqNum, endSeqNum+1)])))
        print (len(AllRepObjs))
        sigReps = printAndParseFiles.SigTEs(AllRepObjs, significantRepeatList)
        print len(sigReps)
        CfamilyMats = FindCandP.combineMatrices(sigReps)
        print ("done with C")
        printAndParseFiles.printListToFile(CfamilyMats, "Calculations_" + organism + "\CFamilyMats.txt")
        PfamilyMats = FindCandP.getPSymm(CfamilyMats)
        printAndParseFiles.printListToFile(PfamilyMats, "Calculations_" + organism + "\PFamilyMats.txt")
        print ("done with P")
    else:
        PfamilyMats = printAndParseFiles.parsePfamilyFile("Calculations_" + organism + "\PFamilyMats.txt")

    # to find global alpha values for R*t
    print ("R")
    RmatDict = FindR.calculateR(PfamilyMats)
    printAndParseFiles.printMatrixToFile(RmatDict, "Calculations_" + organism + "\RFamilyMats.txt" )
       
    # to find global alpha values for m*t
    matSize = 4
    print ("mt")
    mtDict = FindMT.calculateMt(RmatDict, matSize)
    printAndParseFiles.printValuesToFile(mtDict, "Calculations_" + organism + "\MTFamilyVals.txt")
       
    # to find global alpha values for Q
    print ("Q")
    QmatDict = FindQ.calculateQvalues(mtDict, RmatDict, matSize)
     
    printAndParseFiles.printMatrixToFile(QmatDict, "Calculations_" + organism + "\QFamilyMats.txt" )    
    
def RunAlphaGammaTest(chr, organism):
    """
    Runs all Alpha, Gamma test calculations.
    chr: Current chromosome number.
    organism: Organism's name.
    """
    print ("ALPHA GAMMA")
    # to find C alpha, gamma and P alpha, gamma
    seqNames = printAndParseFiles.findSeqNames("geneLocations" + organism + "\GeneLocations%d.txt" % (chr))
    seqCagDictAsym = FindCandP.combineSeqMats("SubMatinGenes" + organism[0].upper() + "\SubMatinGenes%d.txt" % (chr), seqNames)
    seqCagDictSym = copy.deepcopy(seqCagDictAsym)
    printAndParseFiles.printDictionaryListToFile(seqCagDictAsym, "Calculations_" + organism + "\CSeqFamilyMats\CSeqFamilyMats%d.txt" % (chr))
      
    seqPagDictAsym = FindCandP.findPalphaGammaAsymmetric(seqCagDictAsym)
    seqPagDictSym = FindCandP.findPalphaGammaSymmetric(seqCagDictSym)
    printAndParseFiles.printDictionaryListToFile(seqPagDictAsym, "Calculations_" + organism + "\PSeqFamilyMats\PSeqFamilyMatsAsymmetric%d.txt" % (chr))
    printAndParseFiles.printDictionaryListToFile(seqPagDictSym, "Calculations_" + organism + "\PSeqFamilyMats\PSeqFamilyMatsSymmetric%d.txt" % (chr))  
   
    # to find R*t alpha, gamma
    RagDictAsym = FindR.calculateRag(seqPagDictAsym, str(chr), "Asymmetric", "RAssertionFailuresAsymmetric_" + organism + "/")
    RagDictSym = FindR.calculateRag(seqPagDictSym, str(chr), "Symmetric", "RAssertionFailuresSymmetric_" + organism + "/")
    for gene, families in RagDictSym.items():
        for name, Rmatrix in families.items():
            if name not in RagDictAsym[gene]:
                print "From Symmetric"
                print gene, name
                del RagDictSym[gene][name]
    for gene, families in RagDictAsym.items():
        for name, Rmatrix in families.items():
            if name not in RagDictSym[gene]:
                print "From Asymmetric"
                print gene, name
                del RagDictAsym[gene][name]
    printAndParseFiles.printDictionaryMatrixToFile(RagDictAsym, "Calculations_" + organism + "\RSeqFamilyMats\RSeqFamilyMatsAsymmetric%d.txt" % (chr))
    printAndParseFiles.printDictionaryMatrixToFile(RagDictSym, "Calculations_" + organism + "\RSeqFamilyMats\RSeqFamilyMatsSymmetric%d.txt" % (chr))    
  
    # to find Q alpha,gamma and Q gamma
    matSize = 4
    MtagDictAsym = FindMT.calculateMtalphaGamma(RagDictAsym, matSize)
    MtagDictSym = FindMT.calculateMtalphaGamma(RagDictSym, matSize)
    printAndParseFiles.printDictionaryofValuesToFile(MtagDictAsym, "Calculations_" + organism + "\MtSeqFamilyMats\MtSeqFamilyMatsAsymmetric%d.txt" % (chr))
    printAndParseFiles.printDictionaryofValuesToFile(MtagDictSym, "Calculations_" + organism + "\MtSeqFamilyMats\MtSeqFamilyMatsSymmetric%d.txt" % (chr))
    
    QagDictAsym = FindQ.calculateQag(RagDictAsym, MtagDictAsym, matSize)
    QagDictSym = FindQ.calculateQag(RagDictSym, MtagDictSym, matSize)
    printAndParseFiles.printDictionaryMatrixToFile(QagDictAsym, "Calculations_" + organism + "\QSeqFamilyMats\QSeqFamilyMatsAsymmetric%d.txt" % (chr))
    printAndParseFiles.printDictionaryMatrixToFile(QagDictSym, "Calculations_" + organism + "\QSeqFamilyMats\QSeqFamilyMatsSymmetric%d.txt" % (chr))    
    
    QgammaAsym = FindQ.findQgamma(QagDictAsym)
    QgammaSym = FindQ.findQgamma(QagDictSym)
    printAndParseFiles.printMatrixToFile(QgammaAsym, "Calculations_" + organism + "\Qgamma\QgammaAsymmetric%d.txt" % (chr))
    printAndParseFiles.printMatrixToFile(QgammaSym, "Calculations_" + organism + "\Qgamma\QgammaSymmetric%d.txt" % (chr))

    print ("DONE WITH CHR %d" % chr)
    
def RunGammaTest(chr, beginSeqName, globalOrLocalMts, organism):
    """
    Runs all Gamma test calculations.
    chr: Current chromosome number.
    beginSeqName: Beginning of sequence names for query.
    globalOrLocalMts: String identifying whether to use global or local Mt values.
    organism: Organism's name.
    """
    print ("GAMMA")
    # to find P gamma
    familyListAsym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\RSeqFamilyMats\RSeqFamilyMatsAsymmetric%d.txt" % (chr), beginSeqName)#change this to Rag, should be the final area of filtering
    familyListSym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\RSeqFamilyMats\RSeqFamilyMatsSymmetric%d.txt" % (chr), beginSeqName)
    QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "\Qgamma\QgammaAsymmetric%d.txt" % (chr))
    QgammaSym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "\Qgamma\QgammaSymmetric%d.txt" % (chr))
    mtDict = {}
    if (globalOrLocalMts == "Global"):
        mtDict = printAndParseFiles.parseValueFile("Calculations_" + organism + "\MTFamilyVals.txt")
        PalphaGammaHatAsym = FindCandP.calculatePalphaGammaHatGlobal(QgammaAsym, familyListAsym, mtDict)
        PalphaGammaHatSym = FindCandP.calculatePalphaGammaHatGlobal(QgammaSym, familyListSym, mtDict)
        printAndParseFiles.printDictionaryMatrixToFile(PalphaGammaHatAsym, "Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatAsymmetric%d.txt" % (globalOrLocalMts,chr))
        printAndParseFiles.printDictionaryMatrixToFile(PalphaGammaHatSym, "Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatSymmetric%d.txt" % (globalOrLocalMts,chr))

        CalphaGamma = printAndParseFiles.parseCalphaGammaFile("Calculations_" + organism + "\CSeqFamilyMats\CSeqFamilyMats%d.txt" % (chr), beginSeqName)
        
        # calculate this section's L gamma value
        LgammaAsym = FindLg.findLgamma(CalphaGamma, PalphaGammaHatAsym)
        LgammaSym = FindLg.findLgamma(CalphaGamma, PalphaGammaHatSym)
        printAndParseFiles.printValuesToFile(LgammaAsym, "Calculations_" + organism + "\Lgamma%s\LgammaAsymmetric%d.txt" % (globalOrLocalMts,chr))
        printAndParseFiles.printValuesToFile(LgammaSym, "Calculations_" + organism + "\Lgamma%s\LgammaSymmetric%d.txt" % (globalOrLocalMts,chr))
        
    elif (globalOrLocalMts == "Local"):
        mtDictAsym = printAndParseFiles.parseMtalphaGammaFile("Calculations_" + organism + "\MtSeqFamilyMats\MtSeqFamilyMatsAsymmetric%d.txt" % chr, beginSeqName) # create dictionary values for asymmetric and symmetric, perform gamma calculations using them
        PalphaGammaHatAsym = FindCandP.calculatePalphaGammaHatLocal(QgammaAsym, familyListAsym, mtDictAsym)
        mtDictSym = printAndParseFiles.parseMtalphaGammaFile("Calculations_" + organism + "\MtSeqFamilyMats\MtSeqFamilyMatsSymmetric%d.txt" % chr, beginSeqName)
        PalphaGammaHatSym = FindCandP.calculatePalphaGammaHatLocal(QgammaSym, familyListSym, mtDictSym)
        printAndParseFiles.printDictionaryMatrixToFile(PalphaGammaHatAsym, "Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatAsymmetric%d.txt" % (globalOrLocalMts,chr))
        printAndParseFiles.printDictionaryMatrixToFile(PalphaGammaHatSym, "Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatSymmetric%d.txt" % (globalOrLocalMts,chr))

        CalphaGamma = printAndParseFiles.parseCalphaGammaFile("Calculations_" + organism + "\CSeqFamilyMats\CSeqFamilyMats%d.txt" % (chr), beginSeqName)
        
        # calculate this section's L gamma value
        LgammaAsym = FindLg.findLgamma(CalphaGamma, PalphaGammaHatAsym)
        LgammaSym = FindLg.findLgamma(CalphaGamma, PalphaGammaHatSym)
        printAndParseFiles.printValuesToFile(LgammaAsym, "Calculations_" + organism + "\Lgamma%s\LgammaAsymmetric%d.txt" % (globalOrLocalMts,chr))
        printAndParseFiles.printValuesToFile(LgammaSym, "Calculations_" + organism + "\Lgamma%s\LgammaSymmetric%d.txt" % (globalOrLocalMts,chr))

    else:
        print "Error with Global or Local choice. Please type 'Global' or 'Local'"

    print ("DONE WITH CHR %d" % chr)
    
def RunBICTest(chr, beginSeqName, fileName, globalOrLocalTest, organism):
    """
    Runs all BIC test calculations.
    chr: Current chromosome number.
    beginSeqName: Beginning of the sequence name of the queried chromosome. (e.g. "Gene")
    fileName: Output file name.
    globalOrLocalTest: String identifying whether to use global or local Mt values.
    organism: Organism's name.
    """
    print ("BIC")
    if (globalOrLocalTest == "Global"):
        PalphaGammaHatAsym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatAsymmetric%d.txt" % (globalOrLocalTest,chr), beginSeqName)
        PalphaGammaHatSym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatSymmetric%d.txt" % (globalOrLocalTest,chr), beginSeqName)
        LgammaAsym = printAndParseFiles.parseValueFile("Calculations_" + organism + "\Lgamma%s\LgammaAsymmetric%d.txt" % (globalOrLocalTest,chr))
        LgammaSym = printAndParseFiles.parseValueFile("Calculations_" + organism + "\Lgamma%s\LgammaSymmetric%d.txt" % (globalOrLocalTest,chr))
        BICAsym = BIC.findBIC(PalphaGammaHatAsym, LgammaAsym, "A")
        BICSym = BIC.findBIC(PalphaGammaHatSym, LgammaSym, "S")
        printAndParseFiles.compareAndPrintBICs(BICAsym, BICSym, fileName)
    elif (globalOrLocalTest == "Local"):
        PalphaGammaHatAsym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatAsymmetric%d.txt" % (globalOrLocalTest,chr), beginSeqName)
        PalphaGammaHatSym = printAndParseFiles.parsePalphaGamma("Calculations_" + organism + "\PalphaGammaHat%s\PalphaGammaHatSymmetric%d.txt" % (globalOrLocalTest,chr), beginSeqName)
        LgammaAsym = printAndParseFiles.parseValueFile("Calculations_" + organism + "\Lgamma%s\LgammaAsymmetric%d.txt" % (globalOrLocalTest,chr))
        LgammaSym = printAndParseFiles.parseValueFile("Calculations_" + organism + "\Lgamma%s\LgammaSymmetric%d.txt" % (globalOrLocalTest,chr))
        BICAsym = BIC.findBIC(PalphaGammaHatAsym, LgammaAsym, "A")
        BICSym = BIC.findBIC(PalphaGammaHatSym, LgammaSym, "S")
        printAndParseFiles.compareAndPrintBICs(BICAsym, BICSym, fileName)
    else:
        print "Error with Global or Local choice. Please type 'Global' or 'Local'"
        
    print ("DONE WITH CHR %d" % chr)
            
def main():
    """
    File options when executing this file:
    1. Run Global (alpha) test
       - Option number
       - Beginning chromosome
       - Ending chromosome
       - Significant repeat files
       - Organism name
       - .prm file directory (Optional: default is prm_files/)
       - First time to run -> to break up C and P calculations from the rest
    2. Run Local (alpha, gamma) test
       - Option number
       - Beginning chromosome
       - Ending chromosome
       - Organism name
    3. Run Gene Region (gamma) test
       - Option number
       - Beginning chromosome
       - Ending chromosome
       - Beginning of region names
       - Global or local Mts used
       - Organism name
    4. Run BIC test
       - Option number
       - Beginning chromosome
       - Ending chromosome
       - Beginning of region names
       - Global or local Mts used
       - Organism name
    """
    option = int(sys.argv[1])
    beginChr = int(sys.argv[2])
    endChr = int(sys.argv[3])
    if option == 1:
        significantRepeatList = sys.argv[4]
        organism = sys.argv[5]
        if len(sys.argv) == 9:
            prmDirectory = sys.argv[6]
            prmPrefix = sys.argv[7]
            if sys.argv[8] == "True":
                RunGlobalAlphaTest(beginChr, endChr, True, significantRepeatList, organism, prmPrefix, prmDirectory)
            elif sys.argv[8] == "False":
                RunGlobalAlphaTest(beginChr, endChr, False, significantRepeatList, organism, prmPrefix, prmDirectory)
            else:
                print "Error with True or False console parameter. Please type 'True' or 'False'"
        else:
            prmPrefix = sys.argv[6]
            if sys.argv[7] == "True":
                RunGlobalAlphaTest(beginChr, endChr, True, significantRepeatList, organism, prmPrefix)
            elif sys.argv[7] == "False":
                RunGlobalAlphaTest(beginChr, endChr, False, significantRepeatList, organism, prmPrefix)
            else:
                print "Error with True or False console parameter. Please type 'True' or 'False'"
    elif option == 2:
        organism = sys.argv[4]
        for i in range(beginChr, endChr+1):
            RunAlphaGammaTest(i, organism)
    elif option == 3:
        beginSeqName = sys.argv[4] # e.g. "NM_" or "Gene"
        globalOrLocalMts = sys.argv[5]
        organism = sys.argv[6]
        for i in range(beginChr, endChr+1):
            RunGammaTest(i, beginSeqName, globalOrLocalMts) 
    elif option == 4:
        beginSeqName = sys.argv[4] # e.g. "NM_" or "Gene"
        globalOrLocalTest = sys.argv[5]
        organism = sys.argv[6]
        for i in range(beginChr, endChr+1):
            RunBICTest(i, beginSeqName, "Calculations_" + organism + "\BICResults%s\BICResults%d.txt" % (globalOrLocalTest,i), globalOrLocalTest)
            
if __name__== "__main__":
    main()