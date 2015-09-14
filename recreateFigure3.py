'''
recreateFigure3.py
@author: Kristin Helling
'''
import matplotlib.pyplot as plt
import printAndParseFiles
import scipy.stats
import numpy
import sys

def figure3a(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3a in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the C->T substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3a"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#ff0000")
    plt.title("Figure 3a")
    plt.xlabel("$[C \longrightarrow T]_{strand 1} - [C \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()

def figure3b(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3b in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the A->G substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3b"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#ff0000")
    plt.title("Figure 3b")
    plt.xlabel("$[A \longrightarrow G]_{strand 1} - [A \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()
    
def figure3c(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3c in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the G->T substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3c"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#ff0000")
    plt.title("Figure 3c")
    plt.xlabel("$[G \longrightarrow T]_{strand 1} - [G \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()

def figure3d(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3d in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the C->G substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3d"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][1,2]-QgammaAsym[geneNumber][2,1])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#000000")
    plt.title("Figure 3d")
    plt.xlabel("$[C \longrightarrow G]_{strand 1} - [C \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()

def figure3e(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3e in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the A->T substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3e"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][0,3]-QgammaAsym[geneNumber][3,0])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#000000")
    plt.title("Figure 3e")
    plt.xlabel("$[A \longrightarrow T]_{strand 1} - [A \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()

def figure3f(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the histogram shown for Figure 3f in Martin et. al paper.
    Displays the distribution of strand asymmetric regions found  
    for the A->C substitution. This data has not been normalized. 
    Prints out the figure number to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome used (1 by default from main method).
    endChr: Last chromosome used (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure3f"
    histogramPoints = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                histogramPoints.append(QgammaAsym[geneNumber][0,1]-QgammaAsym[geneNumber][3,2])
    n, bins, patches = plt.hist(histogramPoints, bins=30, facecolor="#000000")
    plt.title("Figure 3f")
    plt.xlabel("$[A \longrightarrow C]_{strand 1} - [A \longrightarrow C]_{strand 2}$\n rel. difference")
    plt.ylabel("Distribution")
    plt.show()
    
def main():
    """
    Main method that executes each version of the Figure 3 found in 
    the Martin et. al paper. Redone with both global and local versions.
    The parameters accepted by the console is the final chromosome analyzed
    and the organism name, such as "Maize".
    """
    figure3a("Global", 1, sys.argv[1], sys.argv[2])
    figure3b("Global", 1, sys.argv[1], sys.argv[2])
    figure3c("Global", 1, sys.argv[1], sys.argv[2])
    figure3d("Global", 1, sys.argv[1], sys.argv[2])
    figure3e("Global", 1, sys.argv[1], sys.argv[2])
    figure3f("Global", 1, sys.argv[1], sys.argv[2])
    figure3a("Local", 1, sys.argv[1], sys.argv[2])
    figure3b("Local", 1, sys.argv[1], sys.argv[2])
    figure3c("Local", 1, sys.argv[1], sys.argv[2])
    figure3d("Local", 1, sys.argv[1], sys.argv[2])
    figure3e("Local", 1, sys.argv[1], sys.argv[2])
    figure3f("Local", 1, sys.argv[1], sys.argv[2])
    
if __name__ == '__main__':
    main()