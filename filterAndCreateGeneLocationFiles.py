'''
filterAndCreateGeneLocationFiles.py
@author: Kristin Helling

This file is created to parse through the hg18_gene.mat file to
retrieve the gene locations and capture the number of TE matrices that Martin found.

It is also used to parse through the GFF3 files and have the option to create gene
location files for an organism from these file types.
'''

import sys
import re

def parseMatFile(fp):
    """
    Parses the 'hg18_gene.mat' file to create 'MartinTEsFound.txt',
    which holds the information for each chromosome, its gene number,
    the beginning gene position, the end gene position, and the number
    of TEs found in the gene.
    fp: File pointer for the 'hg18_gene.mat' file.
    """
    countMartinTEs = 0
    geneNumber = 1
    geneBegin = 0
    geneEnd = 0
    chrNum = "chr1"
    fpMartinOut = open("MartinTEsFound.txt", 'w')
    fpMartinOut.write("#ChromsomeNumber\tGeneNumber\tGeneBegin\tGeneEnd\tNumberTEs\n")
    for line in fp:
        line = line.split()
        if geneNumber != int(line[0]) or chrNum != line[1]:
            fpMartinOut.write("%s\t%d\t%d\t%d\t%d\n" % (chrNum,geneNumber,geneBegin,geneEnd,countMartinTEs))
            geneNumber = int(line[0])
            chrNum = line[1]
            countMartinTEs = 0
        if countMartinTEs == 0:
            geneBegin = int(line[2])
            geneEnd = int(line[3])
        countMartinTEs = countMartinTEs + 1
    fpMartinOut.write("%s\t%d\t%d\t%d\t%d\n" % (chrNum,geneNumber,geneBegin,geneEnd,countMartinTEs))
    fpMartinOut.close()
    fp.close()
    
def parseMartinTEFile(fp):
    """
    Parses the 'MartinTEsFound.txt' file to create the
    'GeneLocationsHumanY.txt' files for all of the chromosomes
    in the human genome (where Y is the chromosome number).
    fp: File pointer for the 'MartinTEsFound.txt' file.
    """
    chrNum = 1
    geneNum = 0
    fpChr = open("geneLocations2/GeneLocationsHuman%d.txt" % (chrNum), 'w')
    for line in fp:
        line = line.split()
        if not re.match("^#", line[0]):
            # Filter for 40 or greater TEs found via Martin's findings
            if int(line[4]) >= 40:
                geneNum = geneNum + 1
                if not chrNum == int(line[0].lstrip("chr")):
                    fpChr.close()
                    chrNum = chrNum + 1
                    fpChr = open("geneLocations2/GeneLocationsHuman%d.txt" % (chrNum), 'w')
                    fpChr.write("Gene%d\t%d\t%d\n" % (geneNum,int(line[2]),int(line[3])))
                else:
                    fpChr.write("Gene%d\t%d\t%d\n" % (geneNum,int(line[2]),int(line[3])))  
            elif int(line[4]) < 40:
                print line
    fp.close()
    
def parseGFF3File(fp, organism):
    """
    Parses a GFF3 file to create the 'GeneInformationX.txt'
    file for an organism. (Where X is the organism.)
    fp: File pointer for the GFF3 file.
    organism: The organism name being queried.
    """
    fpOutArr = []
    for line in fp:
        line = line.split()
        if re.match("[1-9]+", line[0]) and line[2] == "gene":
            fpOutArr.append((int(line[0]), int(line[3])-1, line[4], line[6]))  #left closed right open interval, zero based
    fp.close()
    fpOutArr = sorted(fpOutArr, key=lambda gene: gene[0])
    
    count = 1
    fpOut = open("GenomeInformation%s.txt" % organism, 'w')
    fpOut.write("#ChromsomeNumber\tGeneNumber\tGeneBegin\tGeneEnd\tStrand\n")
    for gene in fpOutArr:
        fpOut.write("%s\t%d\t%d\t%s\t%s\n" % ("chr" + str(gene[0]), count, gene[1], gene[2], gene[3]))
        count = count + 1
    fpOut.close()
    
def parseOrganismGeneFile(fp, organism):
    """
    Parses a 'GenomeInformationX.txt' file to create
    all of the 'geneLocationsX/GeneLocationsY.txt' files
    for an organism. (Where X is the organism and Y is the
    chromosome number.)
    fp: File pointer for the 'GenomeInformationX.txt' file.
    organism: The organism name being queried.
    """
    chrNum = 1
    geneNum = 1
    fpChr = open("geneLocations%s/GeneLocations%d.txt" % (organism, chrNum), 'w')
    for line in fp:
        line = line.split()
        if not re.match("^#", line[0]):
            if not chrNum == int(line[0].lstrip("chr")):
                fpChr.close()
                chrNum = chrNum + 1
                fpChr = open("geneLocations%s/GeneLocations%d.txt" % (organism, chrNum), 'w')
                fpChr.write("Gene%s\t%s\t%s\n" % (line[1],line[2],line[3]))
            else:
                fpChr.write("Gene%s\t%s\t%s\n" % (line[1],line[2],line[3]))
            geneNum = geneNum +1
    fp.close()
            
def main():
    """
    File options when executing this file:
    1. Parse through the hg18_gene.mat file to create 'MartinTEsFound.txt'.
    2. Use MartinTEsFound.txt to determine which gene regions were used by Martin et. al.
    3. Parse a GFF3 file to create a 'GenomeInformationX.txt' file based on the queried organism.
    4. Parse the 'GenomeInformationX.txt' file to create all of the gene location files for the organism.
    """
    if sys.argv[1] == "1":
        parseMatFile(open(sys.argv[2]))
    elif sys.argv[1] == "2":
        parseMartinTEFile(open(sys.argv[2]))
    elif sys.argv[1] == "3":
        parseGFF3File(open(sys.argv[2]), sys.argv[3])
    elif sys.argv[1] == "4":
        parseOrganismGeneFile(open(sys.argv[2]), sys.argv[3])
    
if __name__=="__main__":
    main()