# PGine:
# Py/Bioconda software package for calculation of polygenic risk score in plants
# Polygenic risk score
# author: Jana Schwarzerova & Martin Hurta
# pythone v3.10.2

# .............................................................
# Packages:
import pandas as pd

#GWAS
#https://pypi.org/projectpandasgwas/

# .............................................................

#Ukoly
# 1.

# .............................................................
# Loading data:

# Genomic data________________________________________________
gen_info = pd.read_csv(r'Dataset/7c-Genotyping_50K_41722_nG.csv')

# .............................................................
#GWAS

# .............................................................
#Quality Control

# OSETRI!!!

# Sample-level QC: Filter samples with a minimum call rate
sample_call_rate_threshold = 0.95
genotype_data = genotype_data.dropna(thresh=int(sample_call_rate_threshold * len(genotype_data.columns)))

# SNP-level QC: Filter SNPs with a minimum call rate and minor allele frequency
snp_call_rate_threshold = 0.95
maf_threshold = 0.01
genotype_data = genotype_data.dropna(thresh=int(snp_call_rate_threshold * len(genotype_data)), axis=1)
genotype_data = genotype_data.loc[:, (genotype_data.sum(axis=0) / len(genotype_data)) > maf_threshold]

# Additional QC steps can include filtering based on Hardy-Weinberg equilibrium and more

# If necessary, remove duplicate samples
genotype_data = genotype_data.drop_duplicates()

#ROZHODOVACI PRAVIDLO




#implementace LDpred
import numpy as np

# Simulated SNP data (0 = reference allele, 1 = alternate allele)
# Each row represents an individual, and each column represents a SNP
snp_data = np.array([[0, 1, 1, 0, 0],
                     [1, 0, 0, 0, 1],
                     [1, 0, 1, 1, 0],
                     [0, 1, 0, 1, 1]])

# Calculate the allele frequencies
allele_frequencies = snp_data.mean(axis=0)

# Calculate the LD scores for all pairs of SNPs
num_snps = snp_data.shape[1]
ld_scores = np.zeros((num_snps, num_snps))

for i in range(num_snps):
    for j in range(i + 1, num_snps):
        # Calculate LD score (D')
        p_ab = np.mean(snp_data[:, i] * snp_data[:, j])
        d_prime = (p_ab - (allele_frequencies[i] * allele_frequencies[j])) / (
            np.sqrt(allele_frequencies[i] * (1 - allele_frequencies[i]) * allele_frequencies[j] * (1 - allele_frequencies[j]))
        )
        ld_scores[i, j] = ld_scores[j, i] = d_prime

# Print LD scores
print("LD Scores:")
print(ld_scores)



# .............................................................
# PRS




