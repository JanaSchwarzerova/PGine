def genome_wide_association_study(genotype_data, phenotype_data):
    """
    Computes a genome-wide association study, at the end the function returns a p-value
    and effect size alleles

    :param genotype_data: (dataframe) matrix represents SNP information
    :param phenotype_data: (dataframe) represents phenotype information (healthy / sick)
    :return: p-values, effect size (alias effect alleles) and Manhattan plot

    The code was developed using inspiration from PLINK manual and 10.1002/mpr.1608
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
    """

    # Libraries & packages
    import pandas as pd
    import statsmodels.api as sm

    # Load genotype data (SNPs) and phenotype data
    # genotype_data = pd.read_csv(r'FinalTool_GCB_2023\Data\Genotype.csv')
    # phenotype_data = pd.read_csv('FinalTool_GCB_2023\Data\Phenotype.csv')

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
    results_df.to_csv(r"results_GWAS.csv")

    # For signification GWAS:
    #   same GWAS implementation - just with Bonferroni correction

    alpha = 0.05 # Define the significance level (alpha) for Bonferroni correction

    # Calculate the Bonferroni-corrected significance threshold
    num_tests = len(genotype_data_columns)
    bonferroni_threshold = alpha / num_tests

    for snp in genotype_data_columns:
        X = pd.to_numeric(data[snp])
        X = sm.add_constant(X)
        y = data['phenotype']

        model = sm.OLS(y, X).fit()  # Fit linear regression model
        beta = model.params[snp]  # Coefficient represents effect size
        p_value = model.pvalues[snp]  # Calculate p-value

        # Check if the p-value is below the Bonferroni-corrected threshold
        if p_value <= bonferroni_threshold:
            results.append({'SNP': snp, 'P-Value': p_value, 'Effect Size (Beta)': beta})

    significant_results_df  = pd.DataFrame(results)
    significant_results_df .to_csv(r"significant_results_GWAS.csv")


# Other possibility of expansion – GWAS:
# https://pypi.org/projectpandasgwas/


