import os
import random
import csv

# Import modules from the 'initialize' subdirectory
from initialize.pdb_to_fasta import run_fasta_conversion
from initialize.get_matrix_distances import run_distance_matrices
from initialize.get_contact_maps import generate_all_contact_maps


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

    print("\n" + "="*50)
    print(" STARTING PDB PROCESSING PIPELINE")
    print("="*50)

    # STEP 1: Convert PDB to FASTA
    print("STEP 1: Extracting FASTA sequences")
    run_fasta_conversion(PDB_DIR, FASTA_DIR)

    # STEP 2: Compute distance matrices
    print("STEP 2: Computing distance matrices")
    run_distance_matrices(PDB_DIR, DISTANCES_DIR)

    # STEP 3: Generate contact maps
    print("STEP 3: Generating contact maps")
    # Default threshold is 8.0 Å
    generate_all_contact_maps(DISTANCES_DIR, CONTACT_MAPS_DIR)

    print("\n" + "*-"*30 + "*")
    print(" PIPELINE COMPLETED SUCCESSFULLY ")
    print("*-"*30 + "*" + "\n")


if __name__ == "__main__":
    main()