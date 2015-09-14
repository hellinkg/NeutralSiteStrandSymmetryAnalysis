'''
RetrievePlantFamilies.py
@author: Kristin Helling
'''
import re
import glob
import sys

def getZeamaysFamilies(maizeDirectory, maizeFileName):
    """
    Finds all of the families listed within all of the RepeatMasker
    data for Zea mays, based on the information within the .out files.
    maizeDirectory: Complete directory name where the RepeatMasker data is kept
                    for Zea mays.
    maizeFileName: Name used for the prefix of the RepeatMasker data.
    """
    families = []
    for file in glob.glob(maizeDirectory + maizeFileName + "*.fa.out"):
        fp = open(file)
        for line in fp:
            line = line.split()
            if len(line) != 0 and re.match("^[0-9]+",line[0]) and not line[10] in families: #line[10] is repclassfam
                families.append(line[10])
        fp.close()
    printFamilyFile("Zeamays", families)

def getArabidopsisFamilies(arabidopsisDirectory, arabidopsisFileName):
    """
    Finds all of the families listed within all of the RepeatMasker
    data for Arabidopsis thaliana, based on the information within the .out files.
    arabidopsisDirectory: Complete directory name where the RepeatMasker data is kept
                          for Arabidopsis thaliana.
    arabidopsisFileName: Name used for the prefix of the RepeatMasker data.
    """
    families = []
    for file in glob.glob(arabidopsisDirectory + arabidopsisFileName + "*.fa.out"):
        print file
        fp = open(file)
        for line in fp:
            line = line.split()
            if len(line) != 0 and re.match("^[0-9]+",line[0]) and not line[10] in families:
                families.append(line[10])
        fp.close()
    printFamilyFile("Arabidopsisthaliana", families)
        
def printFamilyFile(organism, families):
    """
    Prints out the output from families found within the RepeatMasker data.
    organism: Name of the queried organism.
    families: List of families found from the RepeatMasker data.
    """
    fp = open(organism + "families.txt", "w")
    for i in families:
        fp.write(i + "\n")
    fp.close()
    
def printFamilyInstanceFile(organismDirectory, organism, familyFile):
    """
    Based on the given family file, write out all of the family members
    to a file 'XImportantFamilyInstances.txt', where X is the name of 
    the organism (may be different than the results from the above output).
    organismDirectory: Complete directory name where the RepeatMasker data is kept.
    organism: Name of the queried organism.
    familyFile: 'Xfamilies.txt' file name.
    """
    fp = open(organism + "ImportantFamilyInstances.txt", "w")
    families = [line.rstrip() for line in open(familyFile)]
    instances = []
    for file in glob.glob(organismDirectory + "*.fa.out"):
        print file
        fp2 = open(file)
        for line in fp2:            
            line = line.split()
            if len(line) != 0 and re.match("^[0-9]+",line[0]) and line[10] in families and not line[9] in instances:
                instances.append(line[9])
        fp2.close()
    for i in instances:
        fp.write(i + "\n")
    fp.close()
    
def main():
    """
    File options when executing this file:
    1. Retrieve all of the available families from the Zea mays RepeatMasker data.
    2. Retrieve all of the available families from the Arabidopsis thaliana RepeatMasker data.
    3. Print all of the available family instances from the RepeatMasker data based on the families found from the above operations.
    """
    if sys.argv[1] == "1":
        getZeamaysFamilies(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "2":
        getArabidopsisFamilies(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "3":
        printFamilyInstanceFile(sys.argv[2], sys.argv[3], sys.argv[4])

if __name__ == '__main__':
    main()