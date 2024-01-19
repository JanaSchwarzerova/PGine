# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


#def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
#    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
#if __name__ == '__main__':
#    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/


####################################################################################################
# Libraries & packages
import pandas as pd
import statsmodels.api as sm
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

def genome_wide_association_study_with_cv(genotype_data, phenotype_data):
    """
    Computes a genome-wide association study using 5-fold cross-validation.

    :param genotype_data: (dataframe) matrix represents SNP information
    :param phenotype_data: (dataframe) represents phenotype information (healthy / sick)
    :return: p-values, effect size (alias effect alleles) and Manhattan plot
    """

    # Merge genotype and phenotype data by sample ID
    data = pd.merge(genotype_data, phenotype_data, on='Ind')

    # Initialize lists to store results
    results = []
    significant_results = []
    genotype_data_columns = genotype_data.columns.drop('Ind')

    # Define the number of folds for cross-validation
    num_folds = 5
    kf = KFold(n_splits=num_folds)

    for snp in genotype_data_columns:
        X = pd.to_numeric(data[snp]).values.reshape(-1, 1)
        X_poly = PolynomialFeatures(degree=2).fit_transform(X)  # Use a polynomial of degree 2

        X = sm.add_constant(X_poly)
        y = data['phenotype']

        p_values = []  # Store p-values from each fold

        # Perform 5-fold cross-validation
        for train_index, test_index in kf.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]

            model = LinearRegression().fit(X_train, y_train)  # Fit polynomial regression model on training data
            p_value = model.score(X_test, y_test)  # Calculate R^2 as a measure of goodness-of-fit

            p_values.append(p_value)

        # Calculate the mean R^2 from cross-validation
        mean_p_value = sum(p_values) / num_folds
        beta = model.coef_[1]  # Coefficient for the linear term represents effect size

        results.append({'SNP': snp, 'R^2': mean_p_value, 'Effect Size (Beta)': beta})

        # Check if the mean R^2 is below the Bonferroni-corrected threshold
        alpha = 0.05  # Define the significance level (alpha) for Bonferroni correction
        num_tests = len(genotype_data_columns)
        bonferroni_threshold = alpha / num_tests

        if mean_p_value >= bonferroni_threshold:
            significant_results.append({'SNP': snp, 'R^2': mean_p_value, 'Effect Size (Beta)': beta})

    results_df = pd.DataFrame(results)
    results_df.to_csv("results_GWAS.csv")

    significant_results_df = pd.DataFrame(significant_results)
    significant_results_df.to_csv("significant_results_GWAS.csv")

    # Return both results and significant results dataframes
    return results_df, significant_results_df


# Load genotype data (SNPs) and phenotype data
genotype_data = pd.read_csv(r'C:\Users\xschwa16\Desktop\IGA\Tool\DataBIBM2023\data\Maize_Genotype.csv')
phenotype_data = pd.read_csv(r'C:\Users\xschwa16\Desktop\IGA\Tool\DataBIBM2023\data\Maize_Phenotype.csv')

results, significant_results = genome_wide_association_study_with_cv(genotype_data, phenotype_data)

# .............................................................
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

# doprogramuj vizualizaci


# visualization
plt.boxplot(prs_df)

# Show the plot
plt.tight_layout()
plt.show(block=True)
