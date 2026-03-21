import os

# Importamos los scripts desde la subcarpeta 'initialize'
from initialize.pdb_to_fasta import run_fasta_conversion
from initialize.get_matrix_distances import run_distance_matrices
from initialize.get_contact_maps import generate_all_contact_maps

def main():
    # BASE_DIR es tu carpeta principal 'jazmin'
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    
    # Ruta a la carpeta initialize
    INITIALIZE_DIR = os.path.join(BASE_DIR, "initialize")

    # PDB_DIR ahora apunta correctamente adentro de 'initialize'
    PDB_DIR = os.path.join(INITIALIZE_DIR, "pdb_files")
    
    # Las carpetas de resultados se crearán en tu carpeta principal 'jazmin'
    FASTA_DIR = os.path.join(BASE_DIR, "fasta_files")
    DISTANCES_DIR = os.path.join(BASE_DIR, "real_distances")
    CONTACT_MAPS_DIR = os.path.join(BASE_DIR, "real_contact_maps")

    print("\n" + "="*50)
    print(" INICIANDO PIPELINE DE PROCESAMIENTO PDB")
    print("="*50)

    # PASO 1: PDB a FASTA
    print("\n>>> PASO 1: Extrayendo secuencias FASTA...")
    run_fasta_conversion(PDB_DIR, FASTA_DIR)

    # PASO 2: Matrices de Distancia
    print("\n>>> PASO 2: Calculando matrices de distancia...")
    run_distance_matrices(PDB_DIR, DISTANCES_DIR)

    # PASO 3: Mapas de Contacto
    print("\n>>> PASO 3: Generando mapas de contacto...")
    # El threshold por defecto es 8.0
    generate_all_contact_maps(DISTANCES_DIR, CONTACT_MAPS_DIR)

    print("\n" + "="*50)
    print(" 🎉 PIPELINE COMPLETADO CON ÉXITO 🎉")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()