import numpy as np
import pandas as pd
import statsmodels.api as sm

# Load genotype data (SNPs) and phenotype data
genotype_data = #pd.read_csv('genotype_data.csv')
phenotype_data = #pd.read_csv('phenotype_data.csv')

# Perform quality control on genotype data (e.g., removing low-quality SNPs and samples)
# QC steps might include filtering on missingness, minor allele frequency, Hardy-Weinberg equilibrium, etc.

# Merge genotype and phenotype data by sample ID
data = pd.merge(genotype_data, phenotype_data, on='sample_id')

# Perform association testing for each SNP
results = []
for snp in genotype_data.columns:
    # Fit a linear regression model (or logistic regression for binary traits) for each SNP
    X = data[snp]
    X = sm.add_constant(X)  # Add an intercept term
    y = data['phenotype']

    model = sm.OLS(y, X).fit()  # Replace with logistic regression for binary traits
    p_value = model.pvalues[snp]
    results.append({'SNP': snp, 'P-Value': p_value})
