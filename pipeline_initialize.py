import os
import random
import csv

# Import modules from the 'initialize' subdirectory
from initialize.pdb_to_fasta import run_fasta_conversion
from initialize.get_matrix_distances import run_distance_matrices
from initialize.get_contact_maps import generate_all_contact_maps

# --- DEFINITION OF AMINO ACID CHEMICAL CATEGORIES ---
CATEGORIES = {
    'negatively_charged': ['D', 'E'],
    'positively_charged': ['R', 'K', 'H'],
    'uncharged_polar': ['N', 'Q', 'S', 'T', 'Y'],
    'nonpolar': ['A', 'G', 'V', 'L', 'I', 'P', 'F', 'M', 'W', 'C']
}

# Inverse mapping: amino acid → chemical category
AA_TO_CATEGORY = {aa: cat for cat, aas in CATEGORIES.items() for aa in aas}


def generate_initial_population_by_position(fasta_dir, pop_dir, population_size=200, seed=42):
    """
    Generates an initial population of sequences for a genetic algorithm,
    preserving the chemical category at each position of the original sequence.

    Parameters:
    - fasta_dir: directory containing input FASTA files.
    - pop_dir: directory where the generated population will be stored.
    - population_size: number of sequences to generate per protein.
    - seed: random seed for reproducibility.
    """
    random.seed(seed)
    os.makedirs(pop_dir, exist_ok=True)
    
    # Retrieve all FASTA files in the directory
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(('.fasta', '.fa'))]
    
    if not fasta_files:
        print("  -> No FASTA files found for processing.")
        return

    for file in fasta_files:
        base_name = os.path.splitext(file)[0]
        fasta_path = os.path.join(fasta_dir, file)
        
        # Read original sequence
        original_sequence = []
        with open(fasta_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    original_sequence.append(line.strip().upper())
        
        original_sequence = "".join(original_sequence)
        
        # Generate population
        population = []
        for _ in range(population_size):
            new_sequence = []
            
            # Iterate over each position in the original sequence
            for aa in original_sequence:
                if aa in AA_TO_CATEGORY:
                    # 1. Identify the chemical category at this position
                    current_category = AA_TO_CATEGORY[aa]
                    # 2. Sample a random amino acid within the same category
                    mutated_aa = random.choice(CATEGORIES[current_category])
                    new_sequence.append(mutated_aa)
                else:
                    # Preserve unknown or non-standard amino acids (e.g., 'X')
                    new_sequence.append(aa)
                    
            population.append("".join(new_sequence))
            
        # Save population to CSV
        output_path = os.path.join(pop_dir, f"{base_name}.csv")
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ID', 'Sequence'])
            for i, seq in enumerate(population):
                writer.writerow([f'seq_{i+1}', seq])
        
        print(f"  -> Positional population generated for {base_name} ({population_size} sequences)")


def main():
    # BASE_DIR: root directory of the project
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    
    # Path to the 'initialize' subdirectory
    INITIALIZE_DIR = os.path.join(BASE_DIR, "initialize")

    # Directory containing PDB files
    PDB_DIR = os.path.join(INITIALIZE_DIR, "pdb_files")
    
    # Output directories
    FASTA_DIR = os.path.join(BASE_DIR, "fasta_files")
    DISTANCES_DIR = os.path.join(BASE_DIR, "real_distances")
    CONTACT_MAPS_DIR = os.path.join(BASE_DIR, "real_contact_maps")
    
    # Directory for the genetic algorithm initial population
    POPULATION_DIR = os.path.join(BASE_DIR, "initial_population")

    print("\n" + "="*50)
    print(" STARTING PDB PROCESSING PIPELINE")
    print("="*50)

    # STEP 1: Convert PDB to FASTA
    print("\n>>> STEP 1: Extracting FASTA sequences...")
    run_fasta_conversion(PDB_DIR, FASTA_DIR)

    # STEP 2: Compute distance matrices
    print("\n>>> STEP 2: Computing distance matrices...")
    run_distance_matrices(PDB_DIR, DISTANCES_DIR)

    # STEP 3: Generate contact maps
    print("\n>>> STEP 3: Generating contact maps...")
    # Default threshold is 8.0 Å
    generate_all_contact_maps(DISTANCES_DIR, CONTACT_MAPS_DIR)

    # STEP 4: Generate initial population for the genetic algorithm
    print("\n>>> STEP 4: Generating initial population with positional conservation...")
    generate_initial_population_by_position(
        fasta_dir=FASTA_DIR,
        pop_dir=POPULATION_DIR,
        population_size=200,
        seed=42
    )

    print("\n" + "="*50)
    print(" PIPELINE COMPLETED SUCCESSFULLY ")
    print("="*50 + "\n")


if __name__ == "__main__":
    main()