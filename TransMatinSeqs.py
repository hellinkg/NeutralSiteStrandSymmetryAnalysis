''' 
TransMatinSeqs.py
@author: Kristin Helling 
'''

import sys
from RMfileReader import *
from CountRepeatType import createMatrix
import printAndParseFiles


    
def findRepinLoc(seqIntFile, repObjs, nucNum, outF):
    
    """ 
    Takes the file for the intervals for sequences, the Repeat Objects, and the 
    amount of nucleotides to look outside of the sequence boundaries in order to produce an 
    outfile containing a substitution matrix of each transposable element (TE) found within each designated region.
    seqIntFile: File that contains gene coordinates.
    repObjs: Significant TEs found within the given chromosome of the organism.
    nucNum: The given amount of nucleotides to look previously and following the designated gene boundary.
    outF: The output file name. 
    """
    
    # Gene locations file
    geneLoc = open(seqIntFile)     
    
    # Output: A file containing a substitution matrix for each TE that falls in any of the regions
    repInGene = []
    outFile = open(outF, 'w')
    
    prevCoord = 0
    
    for line in geneLoc:

        lineSp = line.split()
        gene = lineSp[0]
        geneB = int(lineSp[1])
        geneE = int(lineSp[2])

        outFile.write(gene + "\n")
    
        # prepares the gaps for each sequence search
        geneBGap = geneB - nucNum
        geneEGap = geneE + nucNum

        # labels the gene in which one is currently analyzing

        # loop that determines if the TE exists, and how it exists within 
        # the specified region
        te = 0
        while te < len(repObjs):
            outLine = repObjs[te].matchRepeat
            if repObjs[te].posQueryBeg > geneEGap:
                break	    
            elif repObjs[te].posQueryEnd < geneBGap:
                te = te + 1
            elif repObjs[te].posQueryBeg >= geneBGap and repObjs[te].posQueryBeg <= geneEGap:
                repInGene.append(repObjs[te])
                outLine = outLine + " " + createMatrix(repObjs[te], False, False)
                outLine = outLine + "begin TE\n"
                outFile.write(outLine)		
                te = te + 1	    
            elif repObjs[te].posQueryEnd >= geneBGap and repObjs[te].posQueryEnd <= geneEGap:
                repInGene.append(repObjs[te])
                outLine = outLine + " " + createMatrix(repObjs[te], False, False)
                outLine = outLine + "end TE\n"
                outFile.write(outLine)		
                te = te + 1
            elif repObjs[te].posQueryBeg >= geneBGap and repObjs[te].posQueryEnd <= geneEGap:
                repInGene.append(repObjs[te])
                outLine = outLine + " " + createMatrix(repObjs[te], False, False)
                outLine = outLine + "total TE\n"
                outFile.write(outLine)
                te = te + 1
            elif repObjs[te].posQueryBeg < geneBGap and repObjs[te].posQueryEnd > geneEGap:
                repInGene.append(repObjs[te])
                outLine = outLine + " " + createMatrix(repObjs[te], False, False)
                outLine = outLine + "total TE\n"
                outFile.write(outLine)
                te = te + 1

    
    outFile.close()
    geneLoc.close()
    
    
def main():
    """
    Main method for executing TransMatinSeqs.py.
    This program takes in seven arguments in the following order:
    1. Query chromosome prefix name
    2. Gene coordinate file
    3. Chromosome number to query
    4. .prm file directory path
    5. Specific TE file for identification
    6. Nucleotide gap number
    7. Output file name
    """
    # Beginning of chromosome file name (Example: Zea_mays.AGPv3.27.dna.chromosome.)
    chrName = sys.argv[1]
    
    # A file containing a list of intervals on the sequence.  (For example: gene locations.)
    seqIntFile = sys.argv[2]  # "geneLocations/GeneLocationsHuman22.txt"
    
    # RepeatMasker output files for the library created that can parse for the genomic sequence.
    repObjs = unpickleRepeats((chrName + "%d") % int(sys.argv[3]), sys.argv[4])   

    # A file (with full path location included) containing a list of repeat families we are interested in.
    reptFamFile = sys.argv[5]  # "martin_repeats.txt"
    
    # Leeway given for finding TEs around genes
    nucNum = int(sys.argv[6])  # 100
    
    sigAlign = printAndParseFiles.SigTEs(repObjs, reptFamFile)

    outFile = sys.argv[7]
    
    findRepinLoc(seqIntFile, sigAlign, nucNum, outFile)
    
    
if __name__ == "__main__":
    main()
