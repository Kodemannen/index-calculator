import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
from collections import Counter
import pandas as pd
import random
import time
# random.seed(2)

import Levenshtein


from deap import base
from deap import creator
from deap import tools

COLUMNS_OF_INTEREST = ["col_A", "col_B"]

INPUT_FILE = "data.csv"

data_df = pd.read_csv(INPUT_FILE)
dataset_size = len(data_df)
SEQUENCE_LENGTH = len(data_df.iloc[0]["col_A"])


N = 12
popsize = 80            # population size

weigh_median = 2                # for weighing the median more than an individual base
                                # in the fitness function
max_generations = 1000
tournament_size = 4
mutation_rate = 0.2


parents_pool_frac = 1   # create parents_pool_frac*popsize new offspring 
scalefitness = 1        # scaling fitness when used to make probability weights
add_n_randoms = round(0.1*popsize)     # add some randoms each generation

improvement_steps = 100


ignore_last_k_positions = 0


# np.random.seed(0)
# random.seed(0)


# Criterion intervals:
interval = (0.125, 0.4) if N >= 6 else (0, 0.5)       # original:
median_interval = (0.125, 0.4)



def get_base_distribution_for_single_column(df_col, fractional=True):

    # Gets the base distribution for a single column in the dataframe
    # df_col comes in without the header

    distribution = [Counter() for _ in range(SEQUENCE_LENGTH)]

    N_ = len(df_col)


    for sequence in df_col:

        for i, base in enumerate(sequence):
            distribution[i][base] += 1


    for i, position_counter in enumerate(distribution):
        # total = sum(position_counter.values())

        for base in ['A', 'C', 'G', 'T']:
            # .get seems to be a way to avoid key errors by returning 0 if the key is not found
            position_counter[base] = position_counter.get(base, 0) 
            if fractional:
                position_counter[base] /= N_


    distribution_df = pd.DataFrame(distribution).fillna(0)
    distribution_df = distribution_df.reindex(["A", "C", "G", "T"], axis=1)

    return distribution_df


def get_fail_matrices(distribution):

    '''
    fail_matrix is a boolean matrix of the same shape as distribution.values, where
    True means that the value at that position is outside the interval for the specific
    (position, base) combination.
    
    median_fail_matrix is a boolean matrix with four elements, where True means that the median
    distribution for the specific base is outside the interval.
    '''

    fail_matrix = (distribution.values < median_interval[0]) | (distribution.values > median_interval[1])

    medians = np.median(distribution.values, axis=0)
    median_fail_matrix = (medians <= median_interval[0]) | (medians >= median_interval[1])

    return fail_matrix, median_fail_matrix


def criterion_function_for_one_index_column(distribution, n):

    # fail_matrix = (distribution.values <= median_interval[0]) | (distribution.values >= median_interval[1])
    # medians = np.median(distribution.values, axis=0)
    # median_fail_matrix = (medians <= median_interval[0]) | (medians >= median_interval[1])

    fail_matrix, median_fail_matrix = get_fail_matrices(distribution)


    if fail_matrix.any() or median_fail_matrix.any():
        failed = True
    else:
        failed = False

    return failed



def criterion_function():

    pass



def check_if_balanced(sub_df):

    fail_list = []
    for i, header in enumerate(list(sub_df[COLUMNS_OF_INTEREST].columns)): 
        distribution = get_base_distribution_for_single_column(sub_df[header]) 
        failed = criterion_function_for_one_index_column(distribution, len(sub_df))
        print(distribution) 

        fail_list.append(failed)


    return not any(fail_list)






def get_fitnesses(pop, toolbox):
    
    fitnesses = list(map(toolbox.evaluate, pop))


    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    fitnesses = np.array(fitnesses).reshape(-1)
    return fitnesses



def crossover_and_mutate(parent1, parent2, toolbox):
    # Single-point crossover
    point = random.randint(1, 9) # Select a random crossover point

    offspring1 = np.concatenate([parent1[:point], parent2[point:]])
    offspring2 = np.concatenate([parent2[:point], parent1[point:]])

    # Replace duplicates in the offspring with unique indices
    def replace_duplicates(offspring):
        all_indices = set(range(dataset_size))
        used_indices = set(offspring)
        available_indices = list(all_indices - used_indices)
        for i in range(len(offspring)):
            for j in range(i+1, len(offspring)):
                if offspring[i] == offspring[j]:
                    random_index = random.choice(available_indices)
                    offspring[j] = random_index
                    available_indices.remove(random_index)
        return offspring

    offspring1 = replace_duplicates(offspring1)
    offspring2 = replace_duplicates(offspring2)

    # Mutation operation
    def mutate(offspring, mutation_rate=mutation_rate):
        for i in range(len(offspring)):
            if random.random() < mutation_rate:
                all_indices = set(range(dataset_size))
                used_indices = set(offspring)
                available_indices = list(all_indices - used_indices)
                random_index = random.choice(available_indices)
                offspring[i] = random_index
        return offspring

    offspring1 = mutate(offspring1)
    offspring2 = mutate(offspring2)

    offspring1 = creator.Individual(offspring1)
    offspring2 = creator.Individual(offspring2)

    # Evaluate:
    offspring1.fitness.values = toolbox.evaluate(offspring1)
    offspring2.fitness.values = toolbox.evaluate(offspring2)


    return offspring1, offspring2




def get_fitness_from_df(sub_df):

    '''
    Inputs a pandas dataframe with adapter sequences, 
    returns the distance from desired distribution interval
    '''


    distance_metric = 0


    for header in COLUMNS_OF_INTEREST: 

        distribution = get_base_distribution_for_single_column(sub_df[header]) 

        if ignore_last_k_positions:
            distribution = distribution[:-ignore_last_k_positions]


        too_low = distribution.values < interval[0]
        too_high = distribution.values > interval[1]



        distance_too_low  = ((interval[0] - distribution.values)*too_low).sum()
        distance_too_high = ((distribution.values - interval[1])*too_high).sum()


        
        medians = distribution.median(axis=0).values
        medians_too_low = medians < median_interval[0]
        medians_too_high = medians > median_interval[1]


        distance_too_low_medians = ((median_interval[0] - medians)*medians_too_low).sum()
        distance_too_high_medians = ((medians - median_interval[1])*medians_too_high).sum()



        median_distance = distance_too_low_medians + distance_too_high_medians

        total_distance = distance_too_low + distance_too_high + median_distance*weigh_median

        distance_metric += total_distance

    return distance_metric


def evaluate(individual):

    sub_df = data_df.iloc[individual]
    fitness = get_fitness_from_df(sub_df)

    return fitness,



def get_next_generation(pop, 
                        pop_fitnesses, 
                        toolbox):


    current_best_index = np.argmin(pop_fitnesses) 
    current_best = pop[current_best_index]
    current_best_fitness = pop_fitnesses[current_best_index]


    parents_pool_size = round(parents_pool_frac*popsize)

    # round up to nearest even number:
    parents_pool_size = parents_pool_size + parents_pool_size % 2


    fix = 0.0000001
    weights = 1 / (pop_fitnesses**scalefitness + fix)
    weights = weights / weights.sum() # normalize

    parents_indices = np.random.choice(len(pop), size=parents_pool_size, p=weights)

    offspring = []
    offspring_fitnesses = []
    for i in range(0, len(parents_indices), 2):

        parent1 = pop[parents_indices[i]]
        parent2 = pop[parents_indices[i+1]]

        child1, child2 = toolbox.mate(parent1, 
                                      parent2,
                                      toolbox)

        offspring.append(child1)
        offspring.append(child2)

        offspring_fitnesses.append(child1.fitness.values[0])
        offspring_fitnesses.append(child2.fitness.values[0])


    # Add new random individuals to the population:
    randoms = toolbox.population(n=add_n_randoms)
    randoms_fitnesses = get_fitnesses(randoms, toolbox)



    selection_pool = np.concatenate([pop, 
                                     offspring, 
                                     randoms])

    selection_pool_fitnesses = np.concatenate([pop_fitnesses, 
                                               offspring_fitnesses, 
                                               randoms_fitnesses])

    assert len(selection_pool) == len(selection_pool_fitnesses), \
                                    "len(election_pool) != len(selection_pool_fitnesses)"


    weights = 1 / (np.array(selection_pool_fitnesses)**scalefitness + fix)
    weights = weights / weights.sum() # normalize

    new_pop_indices = np.random.choice(len(selection_pool), size=popsize, p=weights)


    # Update population by sampling from selection pool:
    new_pop = np.array(selection_pool)[new_pop_indices]
    new_fitnesses = np.array(selection_pool_fitnesses)[new_pop_indices]



    # replace a random unlucky individual with the current best:
    unlucky_draw = np.random.choice(popsize)
    new_pop[unlucky_draw] = current_best
    new_fitnesses[unlucky_draw] = current_best_fitness

    return new_pop, new_fitnesses


        



def plot_fitness_df():

    fitness_df = pd.read_csv("dev-output/fitness_df_N={}.csv".format(N))
    

    plt.plot(fitness_df["best"], label="best")
    # plt.plot(fitness_df["worst"], label="worst")
    plt.plot(fitness_df["mean"], label="mean")
    # plt.plot(fitness_df["std"], label="std")
    plt.legend()
    plt.title("N={}".format(N))
    # plt.savefig("dev-output/fitness_df_N={}.png".format(N))
    plt.savefig("dev-fig.png".format(N))



def check_solution(N_):

    sub_df = pd.read_csv("dev-output/best_individual_N={}.csv".format(N_))
    balanced = check_if_balanced(sub_df)

    # print(sub_df) 
    if balanced:
        print("Solution is sufficiently balanced!")



def get_distributions(sub_df, fractional=True):
    

    distributions = []
    for i, header in enumerate(COLUMNS_OF_INTEREST): 
        distribution = get_base_distribution_for_single_column(sub_df[header], 
                                                               fractional=fractional) 
        distributions.append(distribution)

    return distributions 



def improve_generation(pop, fitnesses, toolbox, improvement_steps=100):


    for i, individual in enumerate(pop):
        fitness = toolbox.evaluate(individual)[0]


        for j in range(improvement_steps):


            pos_replace = np.random.randint(low=0, high=N)

            new_index = np.random.randint(low=0, high=len(data_df))
            while new_index in individual:
                new_index = np.random.randint(low=0, high=len(data_df))

            new_individual = individual.copy()
            new_individual[pos_replace] = new_index


            new_fitness = toolbox.evaluate(new_individual)[0]


            if new_fitness < fitness:

                fitnesses[i] = new_fitness
                pop[i] = new_individual

                individual = new_individual
                fitness = new_fitness
        

    return pop, fitnesses





def get_levenshtein_distances():


    levenshtein_distances = np.zeros((len(data_df), len(data_df)))
    for header in COLUMNS_OF_INTEREST: 

        # Get the sequences
        seqs = data_df[header].values

        # Get the levenstein distances
        for i in range(len(data_df)):
            for j in range(len(data_df)):

                dist = Levenshtein.distance(seqs[i], seqs[j])

                levenshtein_distances[i, j] += dist


    # Save the levenstein distances

    np.save(f"dev-output/levenstein_distances.npy", levenshtein_distances)





def visualize_levenshtein_dists():

    dists = np.load(f"dev-output/levenstein_distances.npy")

    sns.heatmap(dists)
    plt.show()
    plt.savefig("dev-output/heatmap.pdf")




def build_solution_from_levenshtein_dists():

    dists = np.load(f"dev-output/levenstein_distances.npy")

    '''
    Description:

    Given an MxM matrix that contains the Levenshtein distances from row i to j, do the following:

    0. Start by selecting a random row from the dataset:
    1. Find the row that has the maximum Levenshtein distance from the initial from the first
    2. Find the row that has the maximum Levenshtein distance from the two first
                    :
    n. Find the row that has the maximum Levenshtein distance from the n-1 first

    '''
    

    initial = np.random.randint(low=0, high=dists.shape[0])
    solution_indices = [initial]


    M = len(data_df)

    for i in range(N-1):

        # Choose M random indices, without duplicates, and 
        # without including any of the previously chosen 
        proposals_indices = np.random.choice(dists.shape[0], size=M, replace=False)

        dists_proposals = []


        for j in range(M):

            dist_from_j = 0

            for k in range(len(solution_indices)):
                dist_from_j += dists[solution_indices[k], j] 

            dists_proposals.append(dist_from_j)
                

        farthest_away_index = np.argmax(dists_proposals)
        solution_indices.append(farthest_away_index)

        
    return solution_indices





def main():

    t0 = time.time()
    IND_SIZE = N
    dataset_size = len(data_df)
    
    toolbox = base.Toolbox()

    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    # initialize random individuals:
    toolbox.register("indices", random.sample, range(dataset_size), IND_SIZE)

    # or initialize by choosing rows based on maximum levenshtein distance:
    # toolbox.register("indices", build_solution_from_levenshtein_dists)

    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.indices)


    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    # When registering, tools.initRepeat is called with list and  
    # toolbox.individual as arguments

    # But tools.initRepeat actually takes three arguments (container, func, n),
    # so it's like with functools.partial where it creates a new function that
    # uses container and func as fixed arguments, and then we can input 
    # various values for n



    toolbox.register("evaluate", evaluate)


    toolbox.register("mate", crossover_and_mutate)

    
    pop = toolbox.population(n=popsize)
    

    # Evaluate the initial population
    fitnesses = get_fitnesses(pop, toolbox)
    

    generations = 0


    best_fitness = []
    worst_fitness = []
    mean_fitness = []
    std_fitness = []



    # Begin the evolution
    while generations < max_generations:


        pop, fitnesses =  get_next_generation(pop, 
                                              fitnesses, 
                                              toolbox)

        if improvement_steps > 0:
            pop, fitnesses = improve_generation(pop,
                                                fitnesses,
                                                toolbox,
                                                improvement_steps=improvement_steps)



        best = np.min(fitnesses)
        if best == 0:
            print("Worked")
            break

        best_fitness.append(best)
        worst_fitness.append(np.max(fitnesses))
        mean_fitness.append(np.mean(fitnesses))
        std_fitness.append(np.std(fitnesses))



        
        print("Gen: {}, Best: {:.4f}, Worst: {:.4f}, Mean: {:.4f}, Std: {:.4f}".format(
                                                                             generations, 
                                                                             np.min(fitnesses), 
                                                                             np.max(fitnesses), 
                                                                             np.mean(fitnesses), 
                                                                             np.std(fitnesses)))

        generations += 1


    print("done")
    t1 = time.time() - t0
    print("Time elapsed: {:.2f} seconds".format(t1))

    # Create folder dev-output if not exist:
    if not os.path.exists("dev-output"):
        os.makedirs("dev-output")

    best_individual = data_df.iloc[pop[np.argmin(fitnesses)]]
    best_individual.to_csv(f"dev-output/best_individual_N={N}.csv", index=False)


    
    # create and save a dataframe with the fitness stats
    fitness_df = pd.DataFrame({"best": best_fitness, 
                               "worst": worst_fitness, 
                               "mean": mean_fitness, 
                               "std": std_fitness})

    fitness_df.to_csv(f"dev-output/fitness_df_N={N}.csv", index=False)
        






if __name__=="__main__":

    main()


    # get_levenshtein_distances()
    # visualize_levenshtein_dists()

    # check_solution(N)

    # plot_fitness_df()
    # build_solution_from_levenshtein_dists()


