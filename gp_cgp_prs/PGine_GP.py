"""
PGine: Software for PRS (Polygenic Risk Score) calculation in plants.
Implementation of PRS calculation using Genetic Programming (GP).
gplearn library (https://gplearn.readthedocs.io/en/stable/) is used for GP implementation.

Author: Martin Hurta -  Faculty of Information Technology, Brno University of Technology, Czechia
Version: 1.0.1
Last update: 2024-02-29
"""
# Publicly available libraries and software
import numpy as np
from multiprocessing import Pool, cpu_count
from sklearn.model_selection import KFold
from gplearn.genetic import SymbolicRegressor
from PGine_utils import load_dataset, boxplot_prs, gp_prs


def snp_calculation(data_x, data_y, gp_params):
    """
        Training of GP on data of selected allele.
        5-fold cross validation is used to obtain 5 separate solutions.
        Mean fitness across individual folds is returned together with all 5 GP chromosomes.

    Args:
        data_x (ndarray): 1D array of train input data.
        data_y (ndarray): 1D array of train target data.
        gp_params (dict): Parameters of GP algorithm.

    Returns:
        ndarray: 1D Array with mean fitness over all 5 folds and GP chromosomes from all 5 folds.
    """

    data_x = data_x.reshape(-1, 1)

    # Fitness and found solutions for the current SNP
    fitness = []
    solutions = []

    # Create 5-fold of data for better picture on GP performance
    skf = KFold(5, shuffle=True)

    # Create indices for splitting of data
    splits = skf.split(data_x, data_y)

    # Run GP on each of 5 folds and thus get 5 results and 5 GP solutions
    for i, (train_index, test_index) in enumerate(splits):
        # Get data to be used for training in current fold
        train_x = data_x[train_index]
        train_y = data_y[train_index]

        # Create GP according to provided parameters
        est_gp = SymbolicRegressor(population_size=gp_params["population_size"],
                                   generations=gp_params["generations"],
                                   stopping_criteria=gp_params["stopping_criteria"],
                                   p_crossover=gp_params["p_crossover"],
                                   p_subtree_mutation=gp_params["p_subtree_mutation"],
                                   p_hoist_mutation=gp_params["p_hoist_mutation"],
                                   p_point_mutation=gp_params["p_point_mutation"],
                                   max_samples=gp_params["max_samples"],
                                   verbose=gp_params["verbose"],
                                   parsimony_coefficient=gp_params["parsimony_coefficient"],
                                   random_state=gp_params["random_state"])

        # Run evolution and obtain solution and its fitness on train data
        est_gp.fit(train_x, train_y)

        # Save results of the current fold
        fitness.append(est_gp._program.raw_fitness_)
        solutions.append(est_gp)

    # Return mean train fitness and best solutions over all folds
    return np.append(np.mean(fitness), solutions)


def gp_gwas(genotype_path, phenotype_path, gp_params):
    """
        Learning of GP on all alleles, effectively running GWAS.
        Results of GWAS are saved into results_path. Results include row for each allele and columns:
        allele_name, train fitness and columns for five all parts of CGP genotype.

    Args:
        genotype_path (str): String of path to genotype data.
        phenotype_path (str): String of path to phenotype data.
        gp_params (dict): Parameters of CGP algorithm.
    """

    # Load input SNPs data, target phenotypes and labels of individual SNPs
    data_x, data_y, snp_labels = load_dataset(genotype_path, phenotype_path)

    # Results of GP, axis 0 corresponds to the number of SNPs,
    # axis 1 contains fitness and solution per each fold
    results = np.empty(shape=(len(data_x), 6), dtype=np.object_)

    # Multiprocess results and array of processes
    pool = Pool(processes=cpu_count())
    processes = []

    # Run CGP algorithm on each SNP and get final fitness and solutions
    for snp in data_x:
        processes.append(pool.apply_async(snp_calculation, args=(snp, data_y, gp_params)))

    # Get results of CGP from processes
    for i in range(len(processes)):
        results[i] = processes[i].get()

    return results


def main():

    # Genotype file path
    genotype_path = "/../Data/Genotype.csv"

    # Phenotype file path
    phenotype_path = "/../Data/Phenotype.csv"

    # Parameters of CGP
    gp_params = {
        "population_size": 10,
        "generations": 5,
        "stopping_criteria": 0.01,
        "p_crossover": 0.7,
        "p_subtree_mutation": 0.1,
        "p_hoist_mutation": 0.05,
        "p_point_mutation": 0.1,
        "max_samples": 0.9,
        "verbose": 0,
        "parsimony_coefficient": 0.01,
        "random_state": 0
    }

    # GWAS using CGP
    results = gp_gwas(genotype_path, phenotype_path, gp_params)

    # Calculation of PRS from GWAS data from previous step
    lower_prs, higher_prs = gp_prs(results, genotype_path, phenotype_path)

    # Draw box plots of PRS
    boxplot_prs(lower_prs, higher_prs)


if __name__ == '__main__':
    main()
