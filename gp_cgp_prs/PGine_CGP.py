"""
PGine: Software for PRS (Polygenic Risk Score) calculation in plants.
Implementation of PRS calculation using Cartesian Genetic Programming (CGP).
HAL-CGP library (https://happy-algorithms-league.github.io/hal-cgp/) is used for CGP implementation.

Author: Martin Hurta -  Faculty of Information Technology, Brno University of Technology, Czechia
Version: 1.0.1
Last update: 2024-02-29
"""
import numpy as np
from multiprocessing import Pool, cpu_count
from sklearn.model_selection import KFold
from tqdm import tqdm
import functools
import warnings
import cgp
from PGine_utils import load_dataset, boxplot_prs, gp_prs
import Functions as fc

# Supress warnings created by HAL-CGP library
warnings.filterwarnings("ignore", category=DeprecationWarning)


def objective(individual, train_x, train_y):
    """
        Objective function for evaluation of candidate solutions.
        Fitness of individual is calculated as MSE (Mean Squared Error) on provided data.
        MSE is saved as negative value, because HAL-CGP aims to minimize the fitness.

    Args:
        individual (cgp.individual.IndividualSingleGenome): Individual to calculate the fitness for.
        train_x (numpy.ndarray): Input data for calculation.
        train_y (numpy.ndarray): Labels for input data.

    Returns:
        cgp.individual.IndividualSingleGenome: Individual with updated fitness.

    Raises:
        ValueError: If input data and labels have not the same length.

    """

    if not individual.fitness_is_None():
        return individual

    # Number of inputs and target labels must be the same
    if len(train_x) != len(train_y):
        raise ValueError("Train and test data have different lengths")

    num_of_inputs = len(train_x)  # Number of inputs

    loss = 0.0  # Initial zero loss

    f = individual.to_func()  # Individual candidate solution transformed into function

    # Go through all inputs and labels
    for i in range(num_of_inputs):
        t = train_y[i]
        y = f(train_x[i][0])

        # Calculate squared error and add to overall loss
        loss -= pow(t - y, 2)

    # Calculate mean value of error
    individual.fitness = loss / num_of_inputs

    # Return individual with updated fitness
    return individual


def snp_calculation(data_x, data_y, population_params, genome_params, ea_params, evolve_params):
    """
        Training of CGP on data of selected allele.
        5-fold cross validation is used to obtain 5 separate solutions.
        Mean fitness across individual folds is returned together with all 5 CGP chromosomes.

    Args:
        data_x (numpy.ndarray): 1D array of train input data.
        data_y (numpy.ndarray): 1D array of train target data.
        population_params (dict): Parameters of CGP population.
        genome_params (dict): Parameters of CGP genome.
        ea_params (dict): Parameters of evolutionary algorithm.
        evolve_params (dict): Parameters of evolution process.

    Returns:
        float64: Average fitness obtained on the current SNP.
        numpy.ndarray: CGP chromosomes from all 5 folds.
    """

    k_folds = 5  # Number of folds

    # Reshape data to required format
    data_x = data_x.reshape(-1, 1)

    # Fitness and found solutions for the current SNP
    fitness = np.empty(shape=k_folds, dtype=np.float64)
    solutions = np.empty(shape=k_folds, dtype=np.object_)

    # Create 5-fold of data for better picture on CGP performance
    skf = KFold(k_folds, shuffle=True)

    # Create indices for splitting of data
    splits = skf.split(data_x, data_y)

    # Run CGP on each of 5 folds and thus get 5 results and 5 CGP solutions
    for i, (train_index, test_index) in enumerate(splits):
        # Get data to be used for training in current fold
        train_x = data_x[train_index]
        train_y = data_y[train_index]

        # Create population, evolutionary algorithm and objective functions of CGP by provided parameters
        pop = cgp.Population(**population_params, genome_params=genome_params)
        ea = cgp.ea.MuPlusLambda(**ea_params)
        obj = functools.partial(objective, train_x=train_x, train_y=train_y)

        # Run evolution and obtain solution and its fitness on train data
        cgp.evolve(pop, obj, ea, **evolve_params)

        # Save results of the current fold
        fitness[i] = pop.champion.fitness
        solutions[i] = pop.champion

    # Return mean train fitness and best solutions over all folds
    return np.mean(fitness), solutions


def cgp_gwas(genotype_path, phenotype_path, population_params, genome_params, ea_params, evolve_params):
    """
        Learning of CGP on all alleles, effectively running GWAS.
        Results of GWAS are saved into results_path. Results include row for each allele and columns:
        allele_name, train fitness and columns for five all parts of CGP genotype.

    Args:
        genotype_path (str): String of path to genotype data.
        phenotype_path (str): String of path to phenotype data.
        population_params (dict): Parameters of GP population.
        genome_params (dict): Parameters of GP genome.
        ea_params (dict): Parameters of evolutionary algorithm.
        evolve_params (dict): Parameters of evolution process.
    """

    # Load input SNPs data, target phenotypes and labels of individual SNPs
    data_x, data_y, snp_labels = load_dataset(genotype_path, phenotype_path)

    # Results of CGP, axis 0 corresponds to the number of SNPs,
    # axis 1 contains fitness and solution per each fold
    results = np.empty(shape=(len(data_x), 6), dtype=np.object_)

    # Multiprocess results and array of processes
    pool = Pool(processes=cpu_count())
    processes = []

    # Run CGP algorithm on each SNP and get final fitness and solutions
    for snp in data_x:
        processes.append(pool.apply_async(snp_calculation, args=(snp, data_y, population_params, genome_params,
                                                                 ea_params, evolve_params)))

    # Get results of CGP from processes
    for i in tqdm(range(len(processes))):
        fitness, solutions = processes[i].get()
        results[i][0] = fitness
        results[i][1:] = solutions[:]

    return results


def main():

    # Genotype file path
    genotype_path = "/../Data/Genotype.csv"

    # Phenotype file path
    phenotype_path = "/../Data/Phenotype.csv"

    # Parameters of CGP
    population_params = {"n_parents": 5, "seed": 8188211}

    # Parameters of CGP genomes
    genome_params = {
        "n_inputs": 1,
        "n_outputs": 1,
        "n_columns": 10,
        "n_rows": 5,
        "levels_back": 10,
        "primitives": (cgp.Add, cgp.Sub, cgp.Mul, cgp.ConstantFloat, fc.Identity, fc.Sin, fc.Cos, fc.ProtectedSqrt,
                       fc.Pow, fc.Abs, fc.Minimum, fc.Maximum, fc.ProtectedDiv)
    }

    # Parameters of evolutionary algorithm
    ea_params = {"n_offsprings": 5, "mutation_rate": 0.01, "tournament_size": 2, "n_processes": 1}

    # Parameters of evolution
    evolve_params = {"max_generations": 100, "print_progress": False}

    # GWAS using CGP
    results = cgp_gwas(genotype_path, phenotype_path, population_params, genome_params, ea_params, evolve_params)

    # Calculation of PRS from GWAS data from previous step
    lower_prs, higher_prs = gp_prs(results, genotype_path, phenotype_path)

    # Draw box plots of PRS
    boxplot_prs(lower_prs, higher_prs)


if __name__ == '__main__':
    main()
