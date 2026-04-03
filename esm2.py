#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess
import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import time
from sklearn.metrics import average_precision_score
import threading
from datetime import datetime


# In[ ]:


base_dir = os.getcwd()
esm2_dir = "results_esm2"
results_fasta_dir = os.path.join(esm2_dir, "fasta_files")
genetic_functions_dir = os.path.join(base_dir, "genetic")
initial_pop_dir = os.path.join(base_dir, "initial_population")

log_lock = threading.Lock()
error_log_path = os.path.join(base_dir, "results_esm2", "esm2_errors.log")


# # **Functions**

# Function to create a FASTA file from a sequence

# In[ ]:


def create_fasta(sequence, fasta_results_dir):
    # 1. Use the sequence as the filename
    filename = f"{sequence}.fasta"
    fasta_path = os.path.join(fasta_results_dir, filename)
    
    # 2. Write file using standard FASTA format
    with open(fasta_path, 'w') as f:
        f.write(f">{sequence}\n")  # Header line starts with '>'
        f.write(f"{sequence}\n")   # Sequence in the next line
        
    print(f"FASTA file successfully saved at: {fasta_path}")
    
    # Return the path for downstream usage
    return fasta_path


# Function to predict maps contact with esm2 

# In[ ]:


def generate_contact_map(fasta_path, sequence):
    esm2_dir = "results_esm2"
    contact_maps_dir = os.path.join(esm2_dir, "contact_maps")
    os.makedirs(contact_maps_dir, exist_ok=True)
    
    csv_path = os.path.join(contact_maps_dir, f"{sequence}.csv")
    
    if os.path.exists(csv_path):
        print(f"Skipping: contact map already exists at {csv_path}")
        return csv_path
    
    python_env = "/home/biocomp/anaconda3/envs/esmfold/bin/python"
    script_path = os.path.join("esm2", "get_contact_map.py")
    
    print(f"Running ESM-2 contact map generation for sequence: {sequence}")
    
    try:
        subprocess.run([python_env, script_path, fasta_path, csv_path], check=True)
        print(f"Contact map successfully generated at: {csv_path}")
        return csv_path
        
    except subprocess.CalledProcessError:
        print("Error: ESM-2 script execution failed.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    return None


# Auxiliary functions for fitness

# In[ ]:


def mask_obvious_contacts(csv_path, min_separation=3):
    
    # print(f" Loading matrix and filtering obvious contacts (separation < {min_separation})...")
    
    # 1. Load the CSV using pandas
    df = pd.read_csv(csv_path)
    probability_matrix = df.to_numpy()
    
    # Protein length (L)
    L = probability_matrix.shape[0]
    
    # 2. Create a grid of indices i, j for the entire matrix
    i_indices, j_indices = np.indices((L, L))
    
    # 3. Compute sequence separation: |i - j|
    separation = np.abs(i_indices - j_indices)
    
    # 4. Create the mask: True where separation is very small
    obvious_mask = separation < min_separation

    # 5. Apply the mask: set those probabilities to 0.0
    filtered_matrix = probability_matrix.copy()
    filtered_matrix[obvious_mask] = 0.0
    
    # The lower triangle is also removed to avoid double-counting pairs
    lower_triangle_mask = np.tril(np.ones((L, L)), k=0).astype(bool)
    filtered_matrix[lower_triangle_mask] = 0.0
    
    return filtered_matrix


# # Fitness shannon entropy

# In[ ]:


def fitness_shannon_entropy(filtered_matrix, real_matrix, grid_size=5):

    # 1. Get the protein length (L)
    L = filtered_matrix.shape[0]
    
    if L == 0 or filtered_matrix.shape != real_matrix.shape:
        print("Error: Invalid dimensions or mismatch between matrices.")

    # 2. Get the Top L indices
    flat_predicted = filtered_matrix.flatten()
    
    # np.argsort sorts from smallest to largest, we take the last L and reverse them
    top_l_flat_indices = np.argsort(flat_predicted)[-L:][::-1]
    
    # Convert 1D indices back to 2D coordinates (rows, columns)
    rows, cols = np.unravel_index(top_l_flat_indices, (L, L))
    
    # 3. Identify which of those Top L predictions are correct (True Positives)
    is_true_positive = real_matrix[rows, cols] == 1
    
    # 4. Create the static 5x5 grid
    T_grid = np.zeros((grid_size, grid_size))
    
    # Filter only the coordinates where the model was correct
    tp_rows = rows[is_true_positive]
    tp_cols = cols[is_true_positive]
    
    # 5. Proportional mapping of correct predictions to the 5x5 grid
    for r, c in zip(tp_rows, tp_cols):
        cell_r = (r * grid_size) // L
        cell_c = (c * grid_size) // L
        T_grid[cell_r, cell_c] += 1
        
    # 6. Compute the raw Shannon sum: H = - sum(p_e * ln(p_e))
    H = 0.0
    for count in T_grid.flatten():
        if count > 0:
            # Probability p_e = hits in the cell / L
            p_e = count / L
            H -= p_e * np.log(p_e) 
            
    return H


# # AUC PR fitness

# In[ ]:


def fitness_auc_pr(filtered_matrix, real_matrix, min_separation=3):
    """
    Computes the fitness (AUC_PR) between a predicted probability matrix
    and a real contact map, evaluating only valid contacts.
    """
    L = filtered_matrix.shape[0]
    
    # 1. Extract only valid indices
    # np.triu_indices extracts the upper triangle.
    # The 'k' parameter allows us to skip the main diagonal and nearby neighbors.
    i_indices, j_indices = np.triu_indices(L, k=min_separation)
    
    # 2. Extract the 1D vectors of predictions and true values
    # y_scores: The probabilities from ESM-2 model
    y_scores = filtered_matrix[i_indices, j_indices]
    
    # y_true: The real binary contact map (0s and 1s)
    y_true = real_matrix[i_indices, j_indices]
    
    # Safety check: if there are no true contacts in this range,
    # AUC_PR cannot be computed properly (returns NaN or error).
    if np.sum(y_true) == 0:
        print("Warning: No true contacts in the evaluated range. AUC_PR is undefined.")
    
    # 3. Compute AUC_PR (Average Precision)
    auc_pr = average_precision_score(y_true, y_scores)
    
    return auc_pr


# ## **Genetic algorithm**

# In[ ]:


base_dir = os.getcwd()
genetic_functions_dir = os.path.join(base_dir, "genetic")
initial_pop_dir = os.path.join(base_dir, "initial_population")
import random


# In[ ]:


def tournament(population, k=2):
    contenders = random.sample(population, k)
    winner = max(contenders, key=lambda individual: individual['fitness'])
    return winner


# In[ ]:


POLAR_AA = ['D', 'E', 'R', 'K', 'H', 'N', 'Q', 'S', 'T', 'Y']
NONPOLAR_AA = ['A', 'G', 'V', 'L', 'I', 'P', 'F', 'M', 'W', 'C']

def mutation(sequence, mutation_prob=0.05):
    new_sequence = []
    for aa in sequence:
        if random.random() < mutation_prob:
            if aa in POLAR_AA:
                options = [option for option in POLAR_AA if option != aa]
                new_aa = random.choice(options)
            elif aa in NONPOLAR_AA:
                options = [option for option in NONPOLAR_AA if option != aa]
                new_aa = random.choice(options)
            else:
                new_aa = aa
            new_sequence.append(new_aa)
        else:
            new_sequence.append(aa)
    return "".join(new_sequence)


# In[ ]:


def crossover(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        raise ValueError("Sequences must have the same length.")
    
    length = len(sequence1)
    point1 = random.randint(0, length - 2)
    point2 = random.randint(point1 + 1, length)

    if random.random() < 0.5:
        parent_a, parent_b = sequence1, sequence2
    else:
        parent_a, parent_b = sequence2, sequence1

    child = parent_a[:point1] + parent_b[point1:point2] + parent_a[point2:]

    return child


# In[ ]:


def evaluate_sequence(sequence, target_matrix, fitness_function, cache_dict):

    # 1. In-memory cache: if the sequence was already evaluated, return the stored value
    if sequence in cache_dict:
        # Optional: print("Sequence retrieved from in-memory cache.")
        return cache_dict[sequence]
    
    # 2. Create the FASTA file
    # Ensure the global directory exists before saving
    os.makedirs(results_fasta_dir, exist_ok=True)
    fasta_path = create_fasta(sequence, results_fasta_dir)
    
    # 3. Generate (or retrieve from disk cache) the contact map
    csv_path = generate_contact_map(fasta_path, sequence)
    
    if csv_path is None:
        print(f"Error: Could not generate the contact map for the sequence.")
        fitness_score = -1.0
        cache_dict[sequence] = fitness_score

        with log_lock:
            with open(error_log_path, "a") as log_file:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                mensaje = f"[{timestamp}] ERROR ESM-2 | Secuencia: {sequence} | Motivo: El predictor devolvió None.\n"
                log_file.write(mensaje)

        return fitness_score
        
    # 4. Load and filter the predicted probability matrix (remove diagonals)
    filtered_matrix = mask_obvious_contacts(csv_path)
    
    # 5. Compute the fitness using the function passed as a parameter
    fitness_score = fitness_function(filtered_matrix, target_matrix)
    
    # 6. Store the result in the in-memory cache
    cache_dict[sequence] = fitness_score
    
    return fitness_score


# In[ ]:


'''
import concurrent.futures

def evaluate_batch_parallel(sequences, target_matrix, fitness_function, cache_dict, max_workers=4):
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        futures = {
            executor.submit(evaluate_sequence, seq, target_matrix, fitness_function, cache_dict): seq
            for seq in sequences
        }
        
        for future in concurrent.futures.as_completed(futures):
            seq = futures[future]
            fit = future.result()
            results.append({'sequence': seq, 'fitness': fit})
                
    return results
'''


# In[ ]:


def evaluate_generation_batch(sequences_to_evaluate, target_matrix, fitness_function, cache_dict):
    """
    Takes a list of sequences, creates a multi-FASTA file, calls the ESM-2 environment once,
    and then computes the fitness for each sequence.
    """
    if not sequences_to_evaluate:
        return []

    esm2_dir = "results_esm2"
    contact_maps_dir = os.path.join(esm2_dir, "contact_maps")
    batch_fasta_path = os.path.join(esm2_dir, "current_batch.fasta")
    
    # 1. Create the Multi-FASTA file
    with open(batch_fasta_path, 'w') as f:
        for seq in sequences_to_evaluate:
            f.write(f">{seq}\n")
            f.write(f"{seq}\n")
            
    # 2. Call ESM-2 ONLY once using subprocess
    python_env = "/home/biocomp/anaconda3/envs/esmfold/bin/python"
    script_path = os.path.join("esm2", "get_contact_map_batch.py")
    
    try:
        # Execution will pause here until ESM-2 finishes processing all sequences
        subprocess.run([python_env, script_path, batch_fasta_path, contact_maps_dir], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Fatal error running ESM-2 in batch: {e}")
        # You can handle the error here (e.g., assign negative fitness, etc.)
        
    # 3. Once ESM-2 has saved all CSVs, compute the fitness
    results = []
    for seq in sequences_to_evaluate:
        csv_path = os.path.join(contact_maps_dir, f"{seq}.csv")
        
        if os.path.exists(csv_path):
            filtered_matrix = mask_obvious_contacts(csv_path)
            fitness_score = fitness_function(filtered_matrix, target_matrix)
        else:
            fitness_score = -1.0  # Penalty if prediction failed
            
        cache_dict[seq] = fitness_score
        results.append({'sequence': seq, 'fitness': fitness_score})
        
    return results


# In[ ]:


def genetic_algorithm_esm2(
    target_protein,
    fitness_function,
    generations=30,
    children_count=30,
    p_best=25,
    max_gpu_workers=4,
    alg_folder="resultados"
):

    start_time = time.time()
    fitness_cache = {}
    
    print(f"Starting ESM-2 Genetic Algorithm for: {target_protein}")

    random.seed(26)
    np.random.seed(26)

    # Load target contact map (ground truth)
    print(f"Loading target contact map from real_contact_maps for {target_protein}...")
    target_csv = os.path.join(base_dir, "real_contact_maps", f"{target_protein}.csv")

    if not os.path.exists(target_csv):
        raise FileNotFoundError(f"Target contact map not found at: {target_csv}")

    target_matrix = pd.read_csv(target_csv).to_numpy()

    # Load initial population
    csv_path = os.path.join(initial_pop_dir, f"{target_protein}.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Initial population not found at {csv_path}")

    df = pd.read_csv(csv_path)
    initial_sequences = df['Sequence'].tolist()

    print("Evaluating initial population in parallel (GPU)...")
    population = evaluate_generation_batch(initial_sequences, target_matrix, fitness_function, fitness_cache)

    # Sort and keep unique top individuals
    population = sorted(population, key=lambda x: x['fitness'], reverse=True)
    unique_initial = []
    seen_initial = set()

    for ind in population:
        if ind['sequence'] not in seen_initial:
            unique_initial.append(ind)
            seen_initial.add(ind['sequence'])

    population = unique_initial[:p_best]

    # Initialize fitness cache
    for ind in population:
        fitness_cache[ind['sequence']] = ind['fitness']

    # Store initial statistics
    fitness_vals_init = [ind['fitness'] for ind in population]
    history_stats = [{
        'generation': 0,
        'best_fitness': np.max(fitness_vals_init),
        'mean_fitness': np.mean(fitness_vals_init),
        'median_fitness': np.median(fitness_vals_init)
    }]

    # Evolutionary loop
    for gen in range(generations):

        current_k = 2 if gen < (generations * 2 / 3) else 3
        print(f"\n--- Generation {gen + 1}/{generations} | Tournament size k={current_k} ---")

        # Selection via tournament
        winners_pool = [tournament(population, k=current_k) for _ in range(children_count)]

        # Generate offspring through crossover and mutation
        children_sequences_raw = []
        for _ in range(children_count):
            parent1 = random.choice(winners_pool)
            parent2 = random.choice(winners_pool)

            child_seq = crossover(parent1['sequence'], parent2['sequence'])
            mutated_child = mutation(child_seq)
            children_sequences_raw.append(mutated_child)

        # Identify sequences requiring evaluation
        unique_to_evaluate = [
            seq for seq in children_sequences_raw if seq not in fitness_cache
        ]

        # Retrieve cached fitness values
        children_population = [
            {'sequence': seq, 'fitness': fitness_cache[seq]}
            for seq in children_sequences_raw if seq in fitness_cache
        ]

        # Evaluate new sequences
        if unique_to_evaluate:
            print(f"Sending a batch of {len(set(unique_to_evaluate))} new sequences to the GPU...")
            new_evaluations = evaluate_generation_batch(
                list(set(unique_to_evaluate)), 
                target_matrix, 
                fitness_function, 
                fitness_cache
            )

            children_population.extend(new_evaluations)

            for ind in new_evaluations:
                fitness_cache[ind['sequence']] = ind['fitness']

        # Combine and select next generation
        combined_population = population + children_population
        combined_population = sorted(combined_population, key=lambda x: x['fitness'], reverse=True)

        unique_population = []
        seen_sequences = set()

        for ind in combined_population:
            if ind['sequence'] not in seen_sequences:
                unique_population.append(ind)
                seen_sequences.add(ind['sequence'])

            if len(unique_population) == p_best:
                break

        population = unique_population

        # Compute statistics
        current_fitness_vals = [ind['fitness'] for ind in population]
        best_fit = np.max(current_fitness_vals)
        mean_fit = np.mean(current_fitness_vals)
        median_fit = np.median(current_fitness_vals)

        history_stats.append({
            'generation': gen + 1,
            'best_fitness': best_fit,
            'mean_fitness': mean_fit,
            'median_fitness': median_fit
        })

        print(f"Best: {best_fit:.4f} | Mean: {mean_fit:.4f} | Median: {median_fit:.4f} | Unique: {len(population)}")

    print("\nEvolution completed.")

    end_time = time.time()
    total_time_seconds = end_time - start_time

    # Save results
    esm2_dir = "results_esm2"
    target_dir = os.path.join(esm2_dir, alg_folder)
    os.makedirs(target_dir, exist_ok=True)

    df_results = pd.DataFrame(population)
    save_path = os.path.join(target_dir, f"genetic_results_esm2_{target_protein}.csv")
    df_results.to_csv(save_path, index=False)

    df_stats = pd.DataFrame(history_stats)
    stats_save_path = os.path.join(target_dir, f"statistics_{target_protein}.csv")
    df_stats.to_csv(stats_save_path, index=False)

    time_save_path = os.path.join(target_dir, f"time_{target_protein}.txt")
    with open(time_save_path, "w") as file:
        file.write(f"Sequence: {target_protein}\n")
        file.write(f"Total time: {total_time_seconds:.2f} seconds\n")

    return population


# ## **Execution**

# In[ ]:


target_protein_name = "1s7m"
final_population = genetic_algorithm_esm2(
    target_protein = target_protein_name,
    fitness_function = fitness_auc_pr,      # fitness_auc_pr | fitness_shannon_entropy
    generations = 30,                       # Total number of generations
    children_count = 30,                    # Number of children generated per generation
    p_best = 25,                            # Population size that survives
    max_gpu_workers = 4,                    # Sequences evaluated simultaneously on the GPU
    alg_folder = "hydrophobic_fitness_auc_pr"
)

target_protein_name = "3sb1"
final_population = genetic_algorithm_esm2(
    target_protein = target_protein_name,
    fitness_function = fitness_auc_pr,      # fitness_auc_pr | fitness_shannon_entropy
    generations = 30,                       # Total number of generations
    children_count = 30,                    # Number of children generated per generation
    p_best = 25,                            # Population size that survives
    max_gpu_workers = 4,                    # Sequences evaluated simultaneously on the GPU
    alg_folder = "hydrophobic_fitness_auc_pr"
)

target_protein_name = "1llm" 
final_population = genetic_algorithm_esm2(
    target_protein = target_protein_name,
    fitness_function = fitness_auc_pr,      # fitness_auc_pr | fitness_shannon_entropy
    generations = 30,                       # Total number of generations
    children_count = 30,                    # Number of children generated per generation
    p_best = 25,                            # Population size that survives
    max_gpu_workers = 4,                    # Sequences evaluated simultaneously on the GPU
    alg_folder = "hydrophobic_fitness_auc_pr"
)



