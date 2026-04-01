import pandas as pd
import numpy as np
import os
import glob

def generate_all_contact_maps(input_dir, output_dir, threshold=8.0):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files = glob.glob(os.path.join(input_dir, "*.csv"))
    
    if not files:
        print(f"No files found in {input_dir}")
        return

    print("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    print(f" GENERATING CONTACT MAPS FOR {len(files)} MATRICES")
    print("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")

    for input_path in files:
        file_name = os.path.basename(input_path)
        output_path = os.path.join(output_dir, file_name)
        
        try:
            df_dist = pd.read_csv(input_path)
            df_contact = (df_dist < threshold).astype(int)
            csv_data = df_contact.to_csv(index=False, header=True)
            lines = [line for line in csv_data.split('\n') if line.strip()]
            
            with open(output_path, 'w') as f:
                f.write('\n'.join(lines))

            print(f"Processed map: {file_name}")

        except Exception as e:
            print(f"Error processing {file_name}: {e}")

    print(f"\n Done! All Contact Maps are in {output_dir}/")

if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = os.path.dirname(SCRIPT_DIR)
    generate_all_contact_maps(os.path.join(BASE_DIR, "real_distances"), os.path.join(BASE_DIR, "real_contact_maps"))