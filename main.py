# Py/PGine: software for calculation of polygenic risk score in plants

# author: Jana Schwarzerova
# pythone v3.10.2

# .............................................................
# Libraries & packages
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Load genotype data (SNPs) and phenotype data
genotype_data = pd.read_csv(r'...\Data\Genotype.csv')
phenotype_data = pd.read_csv(r'...\Data\Phenotype.csv')

# .............................................................
# Quality Control
data = genotype_data.T
df = data.drop('Ind', axis = 0)
# ad 3a Checking the amount of SNP errors: Check the data for low-quality SNPs or genotyping errors.
snp_quality_issues = []

for column in df.columns[0:]:  # All columns except 'Marker'
    unique_genotypes = df[column].unique()

    if len(unique_genotypes) > 3:  # We expect only 0, 1, and 2 for genotypes
        snp_quality_issues.append(column)

if snp_quality_issues:
    print("The following markers have SNP errors (low quality or genotyping errors):")
    print(snp_quality_issues)
else:
    print("No markers have SNP errors in the data.")

# ad 3b Allele frequencies: Check that the allele frequencies of SNPs match expectations for your population.
# For simplicity, we will calculate allele frequencies only for reference alleles (0).

allele_frequency_issues = []

for column in df.columns[1:]:  # All columns except 'Marker'
    allele_0_count = (df[column] == 0).sum()
    total_genotypes = len(df[column])
    allele_frequency = allele_0_count / total_genotypes

    if allele_frequency < 0.1 or allele_frequency > 0.9:  # Example: We expect frequencies in the range of 0.1 to 0.9
        allele_frequency_issues.append(column)

if allele_frequency_issues:
    print("The following markers have issues with allele frequencies:")
    print(allele_frequency_issues)
else:
    print("Allele frequencies of SNPs are within the expected range for all markers.")


# .............................................................
# GWAS part
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

# .............................................................
# Calculate the Polygenic Risk Score

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
