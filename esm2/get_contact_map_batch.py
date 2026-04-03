import torch
import esm
import pandas as pd
import sys
import os

def process_batch(multi_fasta_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Read the Multi-FASTA file
    sequences = {}
    current_header = ""
    with open(multi_fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]  # Remove '>'
                sequences[current_header] = ""
            elif current_header:
                sequences[current_header] += line
                
    if not sequences:
        print("No sequences found in the file.")
        return

    # 2. Load the model ONLY ONCE
    print("Loading ESM-2 model (t33_650M)...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    
    if torch.cuda.is_available():
        model = model.cuda()
        
    model.eval()

    # 3. Process each sequence iteratively without reloading the model
    for seq_id, seq in sequences.items():
        csv_out = os.path.join(output_dir, f"{seq_id}.csv")
        
        # Skip if already exists
        if os.path.exists(csv_out):
            continue
            
        print(f"Predicting for: {seq_id}")
        data = [("protein", seq)]
        _, _, batch_tokens = batch_converter(data)
        
        if torch.cuda.is_available():
            batch_tokens = batch_tokens.cuda()
            
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
            
        contacts = results["contacts"][0].cpu().numpy()
        
        # Save CSV
        matrix_size = contacts.shape[0]
        column_names = [f"c{i+1}" for i in range(matrix_size)]
        df_matrix = pd.DataFrame(contacts, columns=column_names)
        df_matrix.to_csv(csv_out, index=False, float_format="%.4f")
        
    print("Batch completed successfully.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python get_contact_map_batch.py <multi_fasta.fasta> <output_dir>")
        sys.exit(1)
        
    multi_fasta_in = sys.argv[1]
    dir_out = sys.argv[2]
    process_batch(multi_fasta_in, dir_out)