"""
PGine: Software for PRS (Polygenic Risk Score) calculation in plants.
Different utilities used by GP and CGP based PRS calculation.

Author: Martin Hurta -  Faculty of Information Technology, Brno University of Technology, Czechia
Version: 1.0.2
Last update: 2024-02-29
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cgp


def load_dataset(genotype_path, phenotype_path, snps=None):
    """
    Load dataset data from provided paths. If list of specific snps is provided, only they will be loaded.

    Args:
        genotype_path (str): String of path to genotype data.
        phenotype_path (str): String of path to phenotype data.
        snps (ndarray): Array of snps names only which will be loaded.
            If none, whole dataset will be loaded.

    Returns:
        ndarray: 2D array of loaded genotype data. Axis 0 corresponds to alleles. Axis 1 corresponds to samples.
        ndarray: 1D array of phenotype data.
        ndarray: 1D array with names of alleles.
    """
    # Load genotype data (SNPs) and phenotype data
    df_genotype = pd.read_csv(genotype_path)
    df_phenotype = pd.read_csv(phenotype_path)

    # If list of snps is provided, filter out the rest
    if snps is not None:
        snps = [snp+1 for snp in snps]
        snps.append(0)
        df_genotype = df_genotype.iloc[:, snps]

    # Connect corresponding data
    df_data = pd.merge(df_genotype, df_phenotype, on='Ind')

    # Remove sample names and create separate input genotype data and phenotypes
    df_data = df_data.drop(['Ind'], axis=1)
    df_data_x = df_data.drop(['phenotype'], axis=1)
    df_data_y = df_data['phenotype']

    # Save SNP labels
    snp_labels = np.asarray(df_data_x.columns)

    # Create numpy arrays for use in calculations
    nd_data_x = np.asarray(df_data_x).transpose()
    nd_data_y = np.asarray(df_data_y)

    return nd_data_x, nd_data_y, snp_labels


def get_the_top_snps(results):
    """
    Get names and functions of SNPs in the best 20% of results obtained by GWAS.

    Args:
        results (ndarray): Array with axis 0 including result for every SNP. Results of individual SNPs contain average
            fitness and five functions (one per each fold).

    Returns:
        ndarray: 1D array with names of SNPs that met the threshold.
        ndarray: 2D array with functions for SNPs that met the threshold.
        ndarray: 1D array with fitness on SNPs that met the threshold.
    """

    fitness = results[:, 0]
    chromosomes = results[:, 1:]

    # Filter the results of SNPs with fitness in the lowest 20% of all SNPs
    n_to_use = int(np.ceil(len(results) / 5))

    p = fitness.argsort()

    fitness = fitness[p]
    chromosomes = chromosomes[p]

    snps = p[:n_to_use]
    filtered_fitness = fitness[:n_to_use]
    filtered_chromosomes = chromosomes[:n_to_use]

    return snps, filtered_chromosomes, filtered_fitness


def boxplot_prs(lower_prs, higher_prs):
    """
    Show boxplot of PRS for samples with phenotype lower or higher than phenotype.

    Args:
        higher_prs (ndarray): PRS values for samples with higher phenotype.
        lower_prs (ndarray): PRS values for samples with lower phenotype.
    """

    # Creation of boxplot
    fig, ax = plt.subplots()

    # Title of plot
    ax.set_title('PRS on data')

    # Creation ox boxplot
    ax.boxplot([lower_prs, higher_prs], labels=['Low', 'High'])

    # Axis labels
    ax.set_xlabel('Rating of plants yield')
    ax.set_ylabel('PRS')

    # Drawing of the boxplot
    plt.show()


def prs_calculation(results, genotype_path, phenotype_path):
    """
    Calculation of PRS from results in provided GWAS file.
        PRS is calculated for all samples with usage of SNPs with fitness in the top 20% of all SNPs.

    Args:
        results (ndarray): Array with axis 0 including result for every SNP. Results of individual SNPs contain
            average fitness and five functions (one per each fold).
        genotype_path (str): String of path to genotype data.
        phenotype_path (str): String of path to phenotype data.

    Raises:
        Exception: Insufficient number of snps with sufficient fitness!
    """

    # Get target phenotype data
    _, nd_data_y, _ = load_dataset(genotype_path, phenotype_path)

    # Calculation of mean phenotype to be used as measure of low vs high plant yield
    mean_phenotype = np.mean(nd_data_y)

    # Creation of masks for higher and lower values
    higher_mask = nd_data_y >= mean_phenotype
    lower_mask = nd_data_y < mean_phenotype

    # Get SNPs in the top 20% and corresponding functions
    snps, chromosomes, fitness = get_the_top_snps(results)

    # Inform if there are not any snps with sufficient threshold
    if len(fitness) < 1:
        raise Exception("Insufficient number of snps with sufficient fitness!")

    # Load dataset, but only selected SNPs
    nd_data_x, nd_data_y, snp_labels = load_dataset(genotype_path, phenotype_path, snps)

    # Results of functions on SNPs data
    results = np.empty(shape=(np.shape(nd_data_x)), dtype=np.object_)

    # Go over each SNP and corresponding function
    for x, snp_chromosomes, i in zip(nd_data_x, chromosomes, range(len(nd_data_x))):

        # Results over all functions
        snp_results = np.empty(shape=(len(snp_chromosomes), len(x)), dtype=np.object_)

        # Calculate SNP by all five chromosomes
        for chromosome, j in zip(snp_chromosomes, range(len(snp_chromosomes))):
            # Using CGP
            if type(chromosome) is cgp.individual.IndividualSingleGenome:
                f = chromosome.to_func()
                res = f(x)
            else:  # Using GP
                res = chromosome.predict(x.reshape(-1, 1))

            snp_results[j] = res

        # Calculate mean over five functions from different folds
        beta_over_snp = np.mean(snp_results, axis=0)

        # Append results over all samples on current SNP
        results[i] = beta_over_snp

    # Calculate sum over SNPs for each sample
    prs_array = np.asarray(np.sum(results, axis=0))

    # Filtering of calculated PRS to ones belonging to samples with phenotype higher or lower than mean
    higher_prs = prs_array[higher_mask]
    lower_prs = prs_array[lower_mask]

    return lower_prs, higher_prs
