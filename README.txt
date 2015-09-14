+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Examining Neutral Substitution Sites to Determine Strand Asymmetry in Organisms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Masters Project, August 2015
Miami University, Oxford, Ohio
Author: Kristin Helling
Advisor: John Karro

This project is used to predict if a region within an organism's genome is 
considered to be strand symmetric or asymmetric based on the respective 
statistical models. It utilizes repeat elements known as transposable 
elements (TEs) to analyze these regions.

*******************************************************************************

Note: Any detailed information about running any of the files are located within
the commented areas of the files. RMfileReader.py contains an argument parser, 
which has a "--h/--help" display for its users.

Process and explanation of these files/the methodology:

1. We ran RepeatMasker, a third-party bioinformatics tool, in order to identify
known repeat elements within our query DNA sequences. (Currently, we have ran 
our tool on the human genome, version hg18, Arabidopsis thaliana, version TAIR10.27,
and Zea mays, version AGPv3.27.)
2. From the RepeatMasker results, we use the .align and .out files to create
Repeat class objects. These class objects are generated via RMfileReader.py, a 
library file created that has been made compatible with Python 2.7 and above. This
library gathers the statistical information and sequence data from both of these 
files to create a custom object that is used for the basis of our analysis. These
resulting objects are stored in .prm (pickled RepeatMasker) files for future use,
via the available pickling technique in Python, which can also be reversed
(unpickled) during our calculations.

Following the previous steps, we have a few other preprocessing steps to 
complete before getting to our main calculations. 

3. In order to perform our local calculations, we need to gather the locations
of the regions in which we will be performing these calculations. These are located
in the "geneLocationsX" directory (described below). To generate these files, use 
filterAndCreateGeneLocationFiles.py. For the human genes, since we were validating
the results from the Martin et. al research, we used the file "hg18_gene.mat" to 
retrieve the same gene locations that they used during their research to generate
the location files. For the other organisms, we retrieved data from the Ensembl
database and parsed the .gff3 files for each of the aforementioned releases of
the plant organisms. This file produces an intermediate file "GenomeInformationX.txt",
where "X" is the organism name.
4. Following this procedure, we produce files called "SubMatinGenesY.txt", where "Y"
is the corresponding chromosome number. These use TransMatinSeqs.py and, for the plant
families, RetrievePlantFamilies.py. The human version already had a generated file for
the designated significant families/family instances that were studied for the organism,
gained from the Martin research called "martin_repeats.txt". RetrievePlantFamilies.py
generates the same information that the "martin_repeats.txt" held, but for either
Arabidopsis thaliana or Zea mays. The formats of these files were basic a basic text file,
with each significant repeat instance on a separate line. TransMatinSeqs.py will then
use the .prm files, significant families/family instances files, and the gene locations
files to generate a transition (or "count") matrix for each significant family instance
found within each region to be analyzed in the queried chromosome. This file will calculate
the transition matrix between the modern sequence and ancestral sequence for any TE found
within the boundaries of the given gene region. Therefore, it has the ability to capture
any TE that begins before the given sequence, but its tail is located within the region,
a TE that begins within the given sequence, but its tail is located outside of the region,
a TE that has its beginning and ending start and end outside of the gene region, but the
region is contained within the TE, as well as a TE that is contained within the region's
boundaries. It will also capture an element and place it in more than one gene region if
it does in fact exist within both gene regions. Also, I have also included the batch files
that were used to automate this process for each organism's genome.

After the preprocessing is completed, we can proceed with our global and local calculations:

Global calculations:
1. Find C_alpha by generating a transition matrix for each TE instance found within the genome,
and sum all of the same family instances into the same family to create a single family
matrix representing all of the instances within the genome.
2. Calculate the probability marix of each family using the symmetric model calculation
to generate the P_alpha matrices for each of these families.
3. Find R_alpha by taking the log of P_alpha. This represents the rate matrix for P_alpha.
4. Calculate Mt_alpha values by taking the sum of the diagonal of the R_alpha matrices, 
and multiplying them by -0.25.
5. Finally, find q_alpha for each family within the genome using the formula R_alpha / Mt_alpha.

Local calculations:
1. Repeat the same calculation steps used for the global values as for the local (alpha, gamma) 
values for both the symmetric and asymmetric models (model calculations found in posted code).
2. Once the above calculations have been repeated, then we calculate q_gamma. This is calculated
using a standard average calculation where each of the q_alpha,gamma for each of the models
is summed, and then divided by the number of families located within the gene region (gamma).
3. Using the q_gamma values, then calculate Phat_alpha,gamma for each model. This is the 
probability matrix based on the distanced value and rate values found within the previous
calculations. These values are calculated using the formula exp(mt*q_gamma), using both local
and global values for mt.
4. Finally calculate the log-likelihood of each model represented by L_gamma. This value 
represents the relationship between the P and C matrices for each family, using the 
summation of each of the elements within the calculated matrices of C_alpha,gamma *
log(Phat_alpha,gamma).
5. After all of these calculations have been performed, calculate the BIC values for
each of the gene regions. Compare these results to one another, and based on the statistical
method for the Bayesian Information Criterion (BIC), the lower of the two values is preferred.
These final results for each chromosome will be printed out to the console for the user that
includes a count of how many genes are predicted to be asymmetric out of the total number of
genes calculated.

*******************************************************************************
Formulas used for our methodology:

Probabilty function:

  P(t) = exp(mqt)

  P(t) -> The discrete Markov chain that represents the probability matrix for substitutions,
          where t is the specific amount of time from the ancestral sequence to the modern 
          sequence. (We do not assume constant or steady-state rate model in our calculations.)
  mt   -> Known as the "distance" to be associated with the count matrix (C).
  q    -> The substitution rate matrix scaled to a single unit.

BIC calculation:

**General form:
  
  BIC = -2 * ln(L) + k * ln(n)

  L -> maximized value of likelihood function
  k -> number of free parameters to be estimated
  n -> number of data points, or the sample size

**Asymmetric model:

  BIC = -2 * L_gamma + (11 + length(Phat_alpha, gamma)) * log(16*length(Phat_alpha,gamma))

**Symmetric model:

  BIC = -2 * L_gamma + (5 + length(Phat_alpha, gamma)) * log(16*length(Phat_alpha,gamma))

*******************************************************************************
Calculation Filters:

1. Select or specify the families in which you would like to study (e.g. "martin_repeats.txt").
2. Each TE family instance must have a transition matrix sum of 40 or more.
3. Within the symmetric model, if either of the lambda values (shown in code) equal zero, filter
   this family from analysis.
4. Each row within the P matrices must be approximately 1.0 (to the .0001).
5. No P matrix should have a value below zero.
6. For each R matrix, the determinant should not be zero.
7. Each R matrix should not have a zero eignvalue.
8. Each R matrix error estimate should be finite and smaller than the designated error tolerance.
9. If any R matrix diagonal value is less than zero, filter this family from analysis.
10. If any R matrix off-diagonal value is greater than -1.0e-5, filter this family from analysis.
11. Each mt value should be finite.
12. Each row within the q matrices must be approximately 0.0 (to the .0001).
13. If any of the rows within a count matrix for alpha, gamma calculations contains zero, filter
    the family from analysis in that gene region.
14. If any R_alpha,gamma matrix exists for one model, but not the other, filter the family from
    analysis in that gene region.
15. If any mt value is equivalent to zero in the Phat calculations, filter the family from analysis
    in that gene region.
16. Phat_alpha,gamma values must not be equal to zero for the L_gamma calculation.
17. L_gamma values will be filtered if they equal zero or NaN.

*******************************************************************************
The hierarchy for the output files for the code is as follows:

mainDirectory
--All .py files
--Calculations_X (Where "X" is the name of the organism studied)
----All alpha (global) values will eventually be printed here
----(after first round of calculations)
----CSeqFamilyMats
------C alpha, gamma matrices
----PSeqFamilyMats
------P alpha, gamma matrices
----PalphaGammaHatGlobal
------P gamma 'hat' matrices, using global Mt values
----PalphaGammaHatLocal
------P gamma 'hat' matrices, using global Mt values
----RSeqFamilyMats
------R alpha, gamma matrices
----MtSeqFamilyMats
------Mt alpha, gamma values
----QSeqFamilyMats
------Q alpha, gamma matrices
----Qgamma
------Q gamma matrices
----LgammaGlobal
------L gamma matrices, using global Mt values
----LgammaLocal
------L gamma matrices, using local Mt values
----BICResultsGlobal
------BIC Results of each gamma, using global Mt values
----BICResultsLocal
------BIC Results of each gamma, using local Mt values
--geneLocationsY (Where "Y" is the first letter of the organism's name)
----All gene locations used for the alpha, gamma values
--SubMatinGenesY (Where "Y" is the first letter of the organism's name)
----All C alpha (local) matrices within each gamma (gene)

Note: Each subdirectory within "Calculations_X" will contain separate files
for both symmetric and asymmetric models.

*******************************************************************************

Addtional files not mentioned in methodology:

Any of the recreateFigureX.py files uploaded to this repository were used to 
generate similar figures as those that were in the Martin et. al paper for
validation purposes for our research.