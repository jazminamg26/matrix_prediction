import os
import glob

# Standard bioinformatics dictionary: 3-letter code -> 1-letter code
AA_DICT = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def convert_pdb_to_fasta(pdb_path, fasta_path):
    sequence = ""
    processed_residues = set()

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM  "):
                chain_id = line[21]
                res_seq = line[22:27].strip()
                res_key = f"{chain_id}_{res_seq}"

                if res_key not in processed_residues:
                    res_name = line[17:20].strip()
                    letter = AA_DICT.get(res_name, 'X')
                    sequence += letter
                    processed_residues.add(res_key)

    base_name = os.path.splitext(os.path.basename(pdb_path))[0]
    with open(fasta_path, 'w') as f:
        f.write(f">{base_name}\n")
        f.write(f"{sequence}\n")

    return len(sequence)

def run_fasta_conversion(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdb_files = glob.glob(os.path.join(input_folder, "*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files found in {input_folder}/")
        return

    print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    print(f" CONVERTING {len(pdb_files)} PDBs TO FASTA")
    print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")

    for pdb_path in pdb_files:
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        fasta_path = os.path.join(output_folder, f"{base_name}.fasta")
        
        try:
            length = convert_pdb_to_fasta(pdb_path, fasta_path)
            print(f"Converted: {base_name}.pdb -> {base_name}.fasta ({length} aa)")
        except Exception as e:
            print(f"Error with {base_name}.pdb: {e}")

    print(f"\n Done! All FASTA files are in {output_folder}/")

if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = os.path.dirname(SCRIPT_DIR)
    run_fasta_conversion(os.path.join(BASE_DIR, "pdb_files"), os.path.join(BASE_DIR, "fasta_files"))