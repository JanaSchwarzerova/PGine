"""
PlottingMaize_ManhattanPlot.py – supporting script for plotting the Manhattan plot
"""
# Libraries & packages
import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt

# Load genotype data (SNPs) and phenotype data
genotype_data = pd.read_csv(r'...\Data\Genotype.csv')
phenotype_data = pd.read_csv(r'...\Data\Phenotype.csv')

# Only if you want to plot the GWAS result using Manhattan Plot:
info_snp = pd.read_csv(r'...\Data\InfoSNP.csv')

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
# Calculate the Polygenic Risk Score
info_snp.rename(columns={'SNP.names': 'SNP'}, inplace=True)
data_manhattan_plot = results_df.merge(info_snp, on='SNP', how='inner')
data_manhattan_plot['NegLogP'] = -np.log10(data_manhattan_plot['P-Value'])

# Load data_manhattan_plot and chromosome boundaries
chromosome_boundaries = [0, 1000000, 2000000]

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

# Initialize x-position and chromosome index
x_position = 0
chromosome_index = 0

# Create the Manhattan plot
for chromosome, group in data_manhattan_plot.groupby('Chromosome'):
    x_values = group['Position'] + x_position
    sns.scatterplot(
        x=x_values,
        y='NegLogP',
        data=group,
        label=f'Chromosome {chromosome}',
        ax=ax,
        s=20,  # Marker size
    )

    # Update x-position for the next chromosome
    x_position += chromosome_boundaries[chromosome_index]
    chromosome_index += 1

# Customize the plot
ax.set_xlabel('Genomic Position')
ax.set_ylabel('-log10(P-Value)')
ax.set_title('Manhattan Plot for GWAS')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()
plt.show()
