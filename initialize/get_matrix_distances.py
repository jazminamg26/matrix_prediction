import os
import glob
import pandas as pd
import numpy as np

def compute_real_distance_matrix(pdb_path, output_path):
    residues = {}
    residue_order = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM  "):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21]
                res_seq = line[22:27].strip()
                
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                res_key = f"{chain_id}_{res_seq}"

                if res_key not in residues:
                    residues[res_key] = {'name': res_name, 'CA': None, 'CB': None}
                    residue_order.append(res_key)

                if atom_name == 'CA':
                    residues[res_key]['CA'] = np.array([x, y, z])
                elif atom_name == 'CB':
                    residues[res_key]['CB'] = np.array([x, y, z])

    representative_coords = []
    for key in residue_order:
        res = residues[key]
        if res['name'] == 'GLY':
            coord = res['CA']
        else:
            coord = res['CB'] if res['CB'] is not None else res['CA']

        if coord is not None:
            representative_coords.append(coord)
        else:
            print(f"Warning: Missing key atoms for residue {key} in {pdb_path}")

    N = len(representative_coords)
    matrix = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            matrix[i, j] = np.linalg.norm(representative_coords[i] - representative_coords[j])

    column_names = [f"c{i+1}" for i in range(N)]
    df_matrix = pd.DataFrame(matrix, columns=column_names)

    csv_text = df_matrix.to_csv(index=False, header=True, float_format="%.4f")
    lines = csv_text.split('\n')
    valid_lines = [line for line in lines if line.strip() != ""]
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(valid_lines))

    return N

def run_distance_matrices(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created folder: {output_folder}/")

    pdb_files = glob.glob(os.path.join(input_folder, "*.pdb"))

    if not pdb_files:
        print(f"No PDB files found in {input_folder}/")
        return

    print("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    print(f" COMPUTING REAL DISTANCES FOR {len(pdb_files)} PDBs")
    print("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")

    for pdb_path in pdb_files:
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        output_path = os.path.join(output_folder, f"{base_name}.csv")

        try:
            length = compute_real_distance_matrix(pdb_path, output_path)
            print(f"✅ Processed: {base_name}.pdb -> {base_name}.csv ({length}x{length} matrix)")
        except Exception as e:
            print(f"Error processing {base_name}.pdb: {e}")

    print(f"\n Done! All CSV files are in {output_folder}/")

if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = os.path.dirname(SCRIPT_DIR)
    run_distance_matrices(os.path.join(BASE_DIR, "pdb_files"), os.path.join(BASE_DIR, "real_distances"))