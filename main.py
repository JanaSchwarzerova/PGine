# PGine: Py/Bioconda software package for calculation of polygenic risk score in plants

# author: Jana Schwarzerova
# pythone v3.10.2

# .............................................................
# Libraries & packages
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Load genotype data (SNPs) and phenotype data
genotype_data = pd.read_csv(r'FinalTool_GCB_2023\Data\Genotype.csv')
phenotype_data = pd.read_csv(r'FinalTool_GCB_2023\Data\Phenotype.csv')

# Only if you want to plot the GWAS result using Manhattan Plot:
info_SNP = pd.read_csv(r'FinalTool_GCB_2023\Data\InfoSNP.csv')

# Merge genotype and phenotype data by sample ID
data = pd.merge(genotype_data, phenotype_data, on='Ind')

# Perform association testing for each SNP
results = []
genotype_data_columns = genotype_data.columns.drop('Ind')

for snp in genotype_data_columns:
    X = pd.to_numeric(data[snp])
    X = sm.add_constant(X)
    y = data['phenotype']
    model = sm.OLS(y, X).fit()  # Fit linear regression model
    beta = model.params[snp]  # Coefficient represents effect size
    p_value = model.pvalues[snp]  # Calculation p-value
    results.append({'SNP': snp, 'P-Value': p_value, 'Effect Size (Beta)': beta})

results_df = pd.DataFrame(results)
effect_sizes = results_df.drop('P-Value', axis=1)

## Calculate the Polygenic Risk Score
genotype_data['Ind'] = effect_sizes['SNP']
genotype_data.rename(columns={'Ind': 'SNP'}, inplace=True)

# Merge genetic_data and effect_sizes DataFrames based on the 'SNP' column
alleles = genotype_data.drop('SNP', axis=1).T

# Calculate the Polygenic Risk Score and store in a list
prs_values = []

for ind in range(0, len(alleles.columns)):
    prs = (alleles.iloc[ind, :] * effect_sizes['Effect Size (Beta)']).sum()
    prs_values.append(prs)

# Create a DataFrame from the prs_values list
prs_df = pd.DataFrame(prs_values)
prs_df.to_csv(r"PolygenicRiskScores.csv")

# Results of PRS
print(prs_df)