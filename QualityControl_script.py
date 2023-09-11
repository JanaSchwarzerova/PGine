"""
QualityControl_script.py – this script serves as an example for simulated data for now
"""

# Quality control (QC) is an important step in the analysis of the Polygenic Risk Score (PRS) from genetic data.
# It helps ensure that your input genetic data is reliable and suitable for PRS analysis.

# The following steps can help you perform a QC analysis as part of a PRS analysis:

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

####################################################################################
# 1. Sample Quality:
# ad 1a: DNA quality control: Check that the DNA has been extracted correctly
#        and is of sufficient concentration and quality.
# ad 1b: Sample quality: Check samples for any signs of contamination or degradation.

# Simulated data with DNA sample concentrations
data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
    'DNA_Concentration': [200, 210, 195, 220, 190]  # For example, in ng/μl
}

# Creating a DataFrame from the data
df = pd.DataFrame(data)

# Plotting DNA concentration for each sample

# Displaying the plot
plt.show()
plt.figure(figsize=(8, 4))
plt.bar(df['Sample'], df['DNA_Concentration'], color='skyblue')
plt.xlabel('Sample')
plt.ylabel('DNA Concentration (ng/μl)')
plt.title('Sample Quality - DNA Concentration')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show(block=True)

# Sample quality check - Comparison with a threshold value
threshold_value = 200  # DNA concentration threshold

# Identifying samples that did not pass the quality check
problematic_samples = df[df['DNA_Concentration'] < threshold_value]

if not problematic_samples.empty:
    print("The following samples did not pass the quality check:")
    print(problematic_samples)
else:
    print("All samples meet quality standards.")


####################################################################################
# 2.Genotypic data quality:
# Checking the amount of genotyping errors:
# ad 2a Check your genotyping data for errors, such as genotyping inconsistencies between samples.
# ad 2b Control of heterozygous calls: Heterozygous calls should be close to the expected 50% for diploid organisms.

# Simulated genotypic data (0, 1, 2 for homozygous, heterozygous, and homozygous for the other allele)
data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
    'SNP1': [0, 1, 2, 2, 1],
    'SNP2': [2, 0, 1, 2, 0],
    # Add more SNPs and samples as needed
}

# Create a DataFrame from the data
df = pd.DataFrame(data)

# 2a. Check for genotyping errors (inconsistencies between samples)
genotyping_errors = []

for column in df.columns[1:]:  # Exclude the 'Sample' column
    unique_genotypes = df[column].unique()

    if len(unique_genotypes) > 1:
        genotyping_errors.append(column)

if genotyping_errors:
    print("The following SNPs have genotyping errors (inconsistencies between samples):")
    print(genotyping_errors)
else:
    print("No genotyping errors found in the data.")

# 2b. Control of heterozygous calls
expected_heterozygosity = 0.5  # Expected proportion of heterozygous calls
heterozygous_issues = []

for index, row in df.iterrows():
    heterozygous_count = sum(1 for genotype in row[1:] if genotype == 1)  # Assuming 1 represents heterozygous
    total_genotypes = len(row) - 1  # Exclude the 'Sample' column
    heterozygous_frequency = heterozygous_count / total_genotypes

    if abs(heterozygous_frequency - expected_heterozygosity) > 0.1:  # Allowing a 10% deviation from the expected 50%
        heterozygous_issues.append(row['Sample'])

if heterozygous_issues:
    print("The following samples have issues with heterozygous calls:")
    print(heterozygous_issues)
else:
    print("Heterozygous calls are within the expected range for all samples.")


####################################################################################
# 3.Quality of markers (SNP Quality):
#
# ad 3a Checking the amount of SNP errors: Check the data for low-quality SNPs or genotyping errors.
# ad 3b Allele frequencies: Check that the allele frequencies of SNPs match expectations for your population.

# Simulated marker data (0, 1, 2 for homozygotes, heterozygotes, and homozygotes for the other allele)
data = {
    'Marker': ['Marker1', 'Marker2', 'Marker3', 'Marker4', 'Marker5'],
    'Sample1': [0, 1, 2, 2, 1],
    'Sample2': [2, 0, 1, 2, 0],
    'Sample3': [1, 1, 2, 0, 0],
    'Sample4': [2, 2, 0, 1, 1],
    'Sample5': [1, 0, 2, 1, 2],
    # Add more markers and samples as needed
}

# Create a DataFrame from the data
df = pd.DataFrame(data)

# ad 3a Checking the amount of SNP errors: Check the data for low-quality SNPs or genotyping errors.
snp_quality_issues = []

for column in df.columns[1:]:  # All columns except 'Marker'
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

####################################################################################
# 4. Correction of population structure:
# If your samples come from different ethnic groups, you should consider correcting for the population structure
# to avoid PRS false positives.

# Simulated genotype data (0, 1, 2 for homozygotes, heterozygotes, and homozygotes for the other allele)
data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
    'SNP1': [0, 1, 2, 2, 1],
    'SNP2': [2, 0, 1, 2, 0],
    'SNP3': [1, 2, 0, 1, 0],
    # Add more SNPs and samples as needed
}

# Create a DataFrame from the data
df = pd.DataFrame(data)

# Extract genotype data (exclude the 'Sample' column)
genotype_data = df.iloc[:, 1:]

# Perform Principal Component Analysis (PCA) for population structure correction
pca = PCA(n_components=2)  # You can adjust the number of components as needed
pca_result = pca.fit_transform(genotype_data)

# Create a new DataFrame for PCA results
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])

# Visualize the PCA results
plt.figure(figsize=(8, 6))
plt.scatter(pca_df['PC1'], pca_df['PC2'])
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA for Population Structure Correction')
plt.grid(True)

# Now it is possible to use the PCA results as covariates in your PRS analysis to correct for population structure.
plt.show(block=True)

####################################################################################
# 5. Data normalization:
# Normalize genotypic data as needed to ensure homogeneity between samples.

# Simulated genotype data (0, 1, 2 for homozygotes, heterozygotes, and homozygotes for the other allele)
data = {
    'Sample': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
    'SNP1': [0, 1, 2, 2, 1],
    'SNP2': [2, 0, 1, 2, 0],
    'SNP3': [1, 2, 0, 1, 0],
    # Add more SNPs and samples as needed
}

# Create a DataFrame from the data
df = pd.DataFrame(data)

# Extract genotype data (exclude the 'Sample' column)
genotype_data = df.iloc[:, 1:]

# Perform standardization to normalize the data
scaler = StandardScaler()
normalized_genotype_data = scaler.fit_transform(genotype_data)

# Create a new DataFrame for the normalized genotype data
normalized_df = pd.DataFrame(data=normalized_genotype_data, columns=genotype_data.columns)
