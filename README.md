# Py/PGine: software for calculation of Polygenic Risk Score (PRS) in plants

> The presented prototype software is implemented in the `main.py` script

## Background
Polygenic Risk Scores have been introduced in human genetics to capture the combined effects of allelic variants on individual outcomes. 
Our objective is to adapt this established procedure to the field of plant biology, thus offering new insights into plant pathophysiology. 
Specifically, this approach opens up possibilities for understanding diseases in individual plants within monoculture farming [1].
Later, we plan to refine and optimize this procedure through the application of Cartesian Genetic Programming (CGP).

## Dataset
The tested dataset was taken from [2]. Thanks to these real data, we simulated a healthy and a sick individual - depending on the grain yield. 

## Our initial prototype consists of two main parts
> 1) Quality Control (QC) data
> 2) PRS calculation

### 1) Quality Control (QC) data
This part of code performs quality control and checks for SNP errors and allele frequencies in genotype data. 
It loads genotype and phenotype data, transposes the genotype data for analysis, and then examines each SNP for genotyping errors, 
flagging SNPs with non-standard genotypes. Additionally, it checks whether the allele frequencies for reference alleles of the SNPs 
fall within the expected range for the population, reporting any issues detected.

### 2) PRS calculation
This second part of code conducts a Genome-Wide Association Study (GWAS) by merging genotype and phenotype data, 
performing association testing for each SNP against the phenotype, and storing the results including p-values and effect sizes. 
It allows calculates Polygenic Risk Scores (PRS) by multiplying allele values by their respective effect sizes and summing them, 
resulting in PRS values for each individual, which are saved in a csv file (as output). 
Finally, it prints and displays the PRS results in a DataFrame.

### Úkol 2




### Reference
[1] Zhu, Y., Chen, H., Fan, J., Wang, Y., Li, Y., Chen, J., Fan, J., Yang, S., Hu, L., Leung, H. and Mew, T.W., 2000. Genetic diversity and disease control in rice. Nature, 406(6797), pp.718-722.

[2] Genome-Wide Analysis of Yield in Europe: Allelic Effects Vary with Drought and Heat Scenarios. 2016. Emilie J. Millet, Claude Welcker, Willem Kruijer, Sandra Negro, Aude Coupel-Ledru, Stéphane D. Nicolas, Jacques Laborde, Cyril Bauland, Sebastien Praud, Nicolas Ranc, Thomas Presterl, Roberto Tuberosa, Zoltan Bedo, Xavier Draye, Björn Usadel, Alain Charcosset, Fred Van Eeuwijk & François Tardieu. 2016. Plant Physiology, 172 (2) 749-764 doi: 10.1104/pp.16.00621
