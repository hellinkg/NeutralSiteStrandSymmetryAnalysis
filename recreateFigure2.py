'''
recreateFigure2.py
@author: Kristin Helling
'''
import matplotlib.pyplot as plt
import printAndParseFiles
import scipy.stats
import numpy
import sys

def getGlobalQGamma(fullFilePath):
    """
    Retrieves the global families and their corresponding values
    to calculate the global Q_gamma value. Returns the Q_gamma
    matrix. (Currently unused.)
    """
    QfamilyMats = printAndParseFiles.parseQgammaFile(fullFilePath) #"Calculations/QFamilyMats.txt"
    globalQgamma = numpy.matrix([[0.0,0.0,0.0,0.0],
                                 [0.0,0.0,0.0,0.0], 
                                 [0.0,0.0,0.0,0.0], 
                                 [0.0,0.0,0.0,0.0]])
    for family, familyMat in QfamilyMats.items():
        globalQgamma += familyMat 
    return numpy.divide(globalQgamma,len(QfamilyMats.items()))
    
def figure2a(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the graph shown for Figure 2a in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    G->T and A->G. Differenciates between gene regions that were
    predicted to be asymmetric versus symmetric gene regions.
    This data has not been normalized. Prints out the figure 
    number and the Pearson correlation coefficient and
    p-value for non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure2a"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                blackPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
                blackPlotPointsY.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1]) 
            else:
                greyPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
                greyPlotPointsY.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',edgecolor="none")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',edgecolor="none")
    ax1.set_xlim([-0.15,0.15])
    ax1.set_ylim([-0.4,0.4])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.07, -0.3), fontsize=20 )
    plt.xlabel("$[G \longrightarrow T]_{strand 1} - [G \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("$[A \longrightarrow G]_{strand 1} - [A \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.title("Figure 2a")
    plt.show()
    print scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))

def figure2b(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the graph shown for Figure 2b in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    G->T and C->T. Differenciates between gene regions that were
    predicted to be asymmetric versus symmetric gene regions.
    This data has not been normalized. Prints out the figure number
    and the Pearson correlation coefficient and p-value for 
    non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure2b"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                blackPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
                blackPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0]) 
            else:
                greyPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
                greyPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',edgecolor="none")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',edgecolor="none")
    ax1.set_xlim([-0.15,0.15])
    ax1.set_ylim([-0.5,0.5])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.07, -0.35), fontsize=20 )
    plt.xlabel("$[G \longrightarrow T]_{strand 1} - [G \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("$[C \longrightarrow T]_{strand 1} - [C \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.title("Figure 2b")
    plt.show()
    print scipy.stats.pearsonr(blackPlotPointsX + greyPlotPointsX, blackPlotPointsY + greyPlotPointsY)
    
def figure2c(globalOrLocal, beginChr, endChr, organism):
    """
    Emulates the graph shown for Figure 2c in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    A->G and C->T. Differenciates between gene regions that were
    predicted to be asymmetric versus symmetric gene regions. 
    This data has not been normalized. Prints out the figure number 
    and the Pearson correlation coefficient and p-value for 
    non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure2c"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    for i in range(beginChr, int(endChr)+1):
        BICResults = printAndParseFiles.parseBICfile("Calculations_" + organism + "/BICResults%s/BICResults%d.txt" % (globalOrLocal,i))
        QgammaAsym = printAndParseFiles.parseQgammaFile("Calculations_" + organism + "/Qgamma/QgammaAsymmetric%d.txt" % i)
        for geneNumber, BICresult in BICResults.items():
            if (BICresult == "Asymmetric"):
                blackPlotPointsX.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
                blackPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0]) 
            else:
                greyPlotPointsX.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
                greyPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',edgecolor="none")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',edgecolor="none")
    ax1.set_xlim([-0.4,0.4])
    ax1.set_ylim([-0.5,0.5])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.2, -0.35), fontsize=20 )
    plt.title("Figure 2c")
    plt.xlabel("$[A \longrightarrow G]_{strand 1} - [A \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.ylabel("$[C \longrightarrow T]_{strand 1} - [C \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.show()
    print scipy.stats.pearsonr(blackPlotPointsX + greyPlotPointsX, blackPlotPointsY + greyPlotPointsY)    
    
def main():
    """
    Main method that executes each version of the Figure 2 found in 
    the Martin et. al paper. Redone with both global and local versions.
    The parameters accepted by the console is the final chromosome analyzed
    and the organism name, such as "Maize".
    """
    #globalQgamma = getGlobalQGamma()
    figure2a("Global", 1, sys.argv[1], sys.argv[2])
    figure2b("Global", 1, sys.argv[1], sys.argv[2])
    figure2c("Global", 1, sys.argv[1], sys.argv[2])
    figure2a("Local", 1, sys.argv[1], sys.argv[2])
    figure2b("Local", 1, sys.argv[1], sys.argv[2])
    figure2c("Local", 1, sys.argv[1], sys.argv[2])
    
if __name__ == '__main__':
    main()