import torch
import esm
import numpy as np
import pandas as pd  # <-- Newly added
import sys
import os

def predict_esm2_contacts(fasta_file, output_csv):
    """
    Predicts a residue–residue contact matrix from a protein sequence
    using the ESM-2 model and saves it as a CSV file.
    """
    
    # 1. Read the FASTA file
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        
    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip()
            
    seq_length = len(sequence)
    
    # 2. Load the ESM-2 model (33-layer version)
    print("Loading ESM-2 model (this may take a few seconds)...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    
    # 3. Move model to GPU if available
    if torch.cuda.is_available():
        model = model.cuda()
        
    model.eval()
    
    # 4. Prepare input data
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()
        
    # 5. Predict contacts using attention-based outputs
    print(f"Predicting contacts for sequence of length {seq_length}...")
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        
    # 6. Extract contact matrix
    contacts = results["contacts"][0].cpu().numpy()
    matrix_size = contacts.shape[0]
    
    # 7. Save matrix with strict formatting (headers, no empty rows)
    column_names = [f"c{i+1}" for i in range(matrix_size)]
    df_matrix = pd.DataFrame(contacts, columns=column_names)
    
    csv_text = df_matrix.to_csv(index=False, header=True, float_format="%.4f")
    
    # === Strict cleanup to remove any empty lines ===
    csv_lines = csv_text.split('\n')
    valid_lines = [line for line in csv_lines if line.strip() != ""]
    final_text = '\n'.join(valid_lines)
    
    with open(output_csv, 'w') as f:
        f.write(final_text)
        
    print(f"Success! ESM-2 matrix ({matrix_size}x{matrix_size}) saved with headers to: {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python generate_esm2_matrix.py <input.fasta> <output.csv>")
        sys.exit(1)
        
    fasta_in = sys.argv[1]
    csv_out = sys.argv[2]
    predict_esm2_contacts(fasta_in, csv_out)