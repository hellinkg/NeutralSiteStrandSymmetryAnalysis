'''
recreateFigure4.py
@author: Kristin Helling
'''
import matplotlib.pyplot as plt
import printAndParseFiles
import scipy.stats
import numpy

def getGeneTranscriptionDirection(beginChr, endChr):
    """
    Using the information found in hg18_gene_strand.bed,
    identify the direction of transcription for each gene analyzed.
    Returns the gene number based on the gene position and its 
    direction of transcription via a dictionary.
    beginChr: First chromosome to be analyzed.
    endChr: Last chromosome to be analyzed (usually the final
            chromosome in the genome).
    """
    fp = open("hg18_gene_strand.bed")
    chrFile = open("geneLocations2/GeneLocationsHuman%d.txt" % beginChr)
    genePositions = {}
    for i in range(beginChr, endChr+1):
        chrFile = open("geneLocations2/GeneLocationsHuman%d.txt" % i)
        genePositions.update( {(line.split()[1],line.split()[2]): (line.split()[0]) for line in chrFile})
        chrFile.close()
    geneTranscriptionDirections = {genePositions[(line.split()[1],line.split()[2])]: line.split()[3] for line in fp if genePositions.has_key((line.split()[1],line.split()[2]))}
    return geneTranscriptionDirections
    
def figure4a(geneTranscriptionDirections, beginChr, endChr):
    """
    Emulates the graph shown for Figure 4a in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    G->T and A->G. Differenciates plot points by identifying 
    transcription directions for each gene region. This data has 
    not been normalized. Prints out the figure number and the 
    Pearson correlation coefficient and p-value for 
    non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure4a"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    QgammaAsym = {}
    for i in range(beginChr, endChr+1):    
        QgammaAsym.update(printAndParseFiles.parseQgammaFile("Calculations/Qgamma/QgammaAsymmetric%d.txt" % i))

    for geneNumber, transcriptionDirection in geneTranscriptionDirections.items():
        if (transcriptionDirection == "+"):
            blackPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
            blackPlotPointsY.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1]) 
        else:
            greyPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
            greyPlotPointsY.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',marker="s")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',marker="+")
    ax1.set_xlim([-0.15,0.15])
    ax1.set_ylim([-0.4,0.4])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.07, -0.3), fontsize=20 )
    plt.xlabel("$[G \longrightarrow T]_{strand 1} - [G \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("$[A \longrightarrow G]_{strand 1} - [A \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.title("Figure 4a")
    plt.show()
    print scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))

def figure4b(geneTranscriptionDirections, beginChr, endChr):
    """
    Emulates the graph shown for Figure 4b in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    G->T and C->T. Differenciates plot points by identifying 
    transcription directions for each gene region. This data has 
    not been normalized. Prints out the figure number and the 
    Pearson correlation coefficient and p-value for 
    non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure4b"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    QgammaAsym = {}
    for i in range(beginChr, endChr+1):    
        QgammaAsym.update(printAndParseFiles.parseQgammaFile("Calculations/Qgamma/QgammaAsymmetric%d.txt" % i))

    for geneNumber, transcriptionDirection in geneTranscriptionDirections.items():
        if (transcriptionDirection == "+"):
            blackPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
            blackPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0]) 
        else:
            greyPlotPointsX.append(QgammaAsym[geneNumber][2,3]-QgammaAsym[geneNumber][1,0])
            greyPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',marker="s")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',marker="+")
    ax1.set_xlim([-0.15,0.15])
    ax1.set_ylim([-0.5,0.5])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.07, -0.35), fontsize=20 )
    plt.xlabel("$[G \longrightarrow T]_{strand 1} - [G \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.ylabel("$[C \longrightarrow T]_{strand 1} - [C \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.title("Figure 4b")
    plt.show()
    print scipy.stats.pearsonr(blackPlotPointsX + greyPlotPointsX, blackPlotPointsY + greyPlotPointsY)
    
def figure4c(geneTranscriptionDirections, beginChr, endChr):
    """
    Emulates the graph shown for Figure 4c in Martin et. al paper.
    Displays the correlation between strand asymmetry rates for 
    A->G and C->T. Differenciates plot points by identifying 
    transcription directions for each gene region. This data has 
    not been normalized. Prints out the figure number and the 
    Pearson correlation coefficient and p-value for 
    non-correlation to the console.
    globalOrLocal: Determines the distance version used (mt_alpha or 
                   mt_alpha,gamma) for these calculations.
    beginChr: First chromosome analyzed (1 by default from main method).
    endChr: Last chromosome analyzed (usually last chromosome in organism).
    organism: Organism rates being graphed/analyzed.
    """
    print "figure4c"
    blackPlotPointsX = []
    blackPlotPointsY = []
    greyPlotPointsX = []
    greyPlotPointsY = []
    QgammaAsym = {}
    for i in range(beginChr, endChr+1):
        QgammaAsym.update(printAndParseFiles.parseQgammaFile("Calculations/Qgamma/QgammaAsymmetric%d.txt" % i))

    for geneNumber, transcriptionDirection in geneTranscriptionDirections.items():
        if (transcriptionDirection == "+"):
            blackPlotPointsX.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
            blackPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0]) 
        else:
            greyPlotPointsX.append(QgammaAsym[geneNumber][0,2]-QgammaAsym[geneNumber][3,1])
            greyPlotPointsY.append(QgammaAsym[geneNumber][1,3]-QgammaAsym[geneNumber][2,0])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(greyPlotPointsX, greyPlotPointsY,c='#b0b0b0',marker="s")
    ax1.scatter(blackPlotPointsX, blackPlotPointsY,c='#000000',marker="+")
    ax1.set_xlim([-0.4,0.4])
    ax1.set_ylim([-0.5,0.5])
    ax1.set_xticks([0.0], minor=True)
    ax1.set_yticks([0.0], minor=True)
    ax1.yaxis.grid(True, which='minor')
    ax1.xaxis.grid(True, which='minor')
    ax1.annotate("r =" + str(round(scipy.stats.pearsonr((blackPlotPointsX + greyPlotPointsX), (blackPlotPointsY + greyPlotPointsY))[0],2)), (0.2, -0.35), fontsize=20 )
    plt.title("Figure 4c")
    plt.xlabel("$[A \longrightarrow G]_{strand 1} - [A \longrightarrow G]_{strand 2}$\n rel. difference")
    plt.ylabel("$[C \longrightarrow T]_{strand 1} - [C \longrightarrow T]_{strand 2}$\n rel. difference")
    plt.show()
    print scipy.stats.pearsonr(blackPlotPointsX + greyPlotPointsX, blackPlotPointsY + greyPlotPointsY)
    
def main():
    """
    Main method that executes each version of the Figure 4 found in 
    the Martin et. al paper. Redone with both global and local versions.
    This figure has only been evaluated based on the human data for this
    project because of the low correlation values currently found for the
    plant analysis performed (r <= 0.02).
    """
    geneTranscriptionDirection = getGeneTranscriptionDirection(1,22)
    figure4a(geneTranscriptionDirection, 1, 22)
    figure4b(geneTranscriptionDirection, 1, 22)
    figure4c(geneTranscriptionDirection, 1, 22)
    
if __name__ == '__main__':
    main()