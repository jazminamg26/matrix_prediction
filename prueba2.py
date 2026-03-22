import os
from Bio import SeqIO

def get_fasta_lengths():
    # Buscamos todos los archivos que terminen en .fasta en la carpeta actual
    fasta_files = [f for f in os.listdir('.') if f.endswith('.fasta')]
    
    if not fasta_files:
        print("No se encontraron archivos .fasta en este directorio.")
        return

    print(f"{'Archivo':<15} | {'Longitud (aa)':<15}")
    print("-" * 35)

    for file in sorted(fasta_files):
        try:
            for record in SeqIO.parse(file, "fasta"):
                # Imprime el nombre del archivo y la longitud de la secuencia
                print(f"{file:<15} | {len(record.seq):<15}")
        except Exception as e:
            print(f"Error procesando {file}: {e}")

if __name__ == "__main__":
    get_fasta_lengths()