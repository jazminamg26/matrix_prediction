import subprocess
import os
import argparse
import sys

# --- RUTAS INTELIGENTES (Dinámicas) ---
BASE_DIR = os.path.abspath(os.getcwd()) 
TRROSETTA_DIR = os.path.join(BASE_DIR, "trRosetta") 

# 1. Añadimos trRosetta al "path" para poder importar tu otro script
sys.path.append(TRROSETTA_DIR)

try:
    from generate_msa_with_colabfold import generar_msa_y_limpiar
except ImportError:
    print("❌ Error: No se pudo importar 'generate_msa_with_colabfold.py'.")
    print(f"Asegúrate de que el archivo esté dentro de la carpeta: {TRROSETTA_DIR}")
    sys.exit(1)

def sanitizar_a3m(file_path):
    """Filtro purificador: Limpia el A3M para evitar que trRosetta colapse"""
    if not os.path.exists(file_path):
        return
    
    with open(file_path, 'r') as f:
        lineas = f.readlines()
        
    lineas_limpias = []
    longitud_referencia = 0
    
    for linea in lineas:
        linea = linea.strip()
        # Ignorar líneas vacías y COMENTARIOS de ColabFold
        if not linea or linea.startswith("#"):
            continue 
            
        if linea.startswith(">"):
            lineas_limpias.append(linea)
        else:
            seq = ''.join([c for c in linea if not c.islower()])
            
            if longitud_referencia == 0:
                longitud_referencia = len(seq)
                
            if len(seq) == longitud_referencia:
                lineas_limpias.append(seq)
                
    with open(file_path, 'w') as f:
        f.write('\n'.join(lineas_limpias) + '\n')

def ejecutar_pipeline_trRosetta(fasta_filename, skip_msa=False):
    nombre_base = os.path.splitext(os.path.basename(fasta_filename))[0]
    
    # Rutas de entrada y salida
    fasta_path = os.path.join(BASE_DIR, "fasta_files", fasta_filename)
    
    carpeta_msa_base = os.path.join(TRROSETTA_DIR, "msa")
    msa_dir = os.path.join(carpeta_msa_base, nombre_base)
    a3m_file = os.path.join(msa_dir, f"{nombre_base}.a3m")
    
    dir_resultados = os.path.join(BASE_DIR, "resultados_trRosetta")
    npz_file = os.path.join(dir_resultados, f"{nombre_base}_prediccion.npz")
    csv_file = os.path.join(dir_resultados, f"{nombre_base}_distancias.csv")

    # Regresamos a la ruta específica del entorno que tiene TensorFlow 1.x
    PYTHON_TRROSETTA = "/home/biocomp/anaconda3/envs/trrosetta_env/bin/python"
    
    # Nombre de script actualizado a lo que solicitaste
    SCRIPT_EXTRAER = os.path.join(TRROSETTA_DIR, "distance_matrix.trRosetta.py")

    os.makedirs(dir_resultados, exist_ok=True)

    print("=================================================")
    print(f" 🧬 INICIANDO PIPELINE trRosetta PARA: {nombre_base}")
    print("=================================================")

    # === LÓGICA INTELIGENTE: AUTO-DETECTAR ALINEAMIENTO ===
    if os.path.exists(a3m_file) and not skip_msa:
        print(f"\n💡 ¡Alineamiento previo detectado! ({a3m_file})")
        print("Activando automáticamente el salto de ColabFold para ahorrar tiempo...")
        skip_msa = True

    # === PASOS 1 Y 2: COLABFOLD Y LIMPIEZA ===
    if not skip_msa:
        print("\n[1/4 y 2/4] Generando MSA y limpiando (usando módulo externo)...")
        if not os.path.exists(fasta_path):
            print(f"❌ Error: No se encuentra {fasta_path}")
            sys.exit(1)
            
        generar_msa_y_limpiar(fasta_path, carpeta_msa_base)
        
    else:
        print(f"\n⏭️ [1/4 y 2/4] Saltando ColabFold. Usando alineamiento previo: {a3m_file}")
        if not os.path.exists(a3m_file):
            print(f"❌ Error: No se encontró el archivo {a3m_file}")
            sys.exit(1)

    print("      Purificando el archivo .a3m para compatibilidad estricta...")
    sanitizar_a3m(a3m_file)
    print(f"✅ Archivo {a3m_file} listo, limpio y sanitizado.")

    # === PASO 3: PREDICCIÓN DE RED NEURONAL ===
    print("\n[3/4] Prediciendo distancias (Red Neuronal)...")
    try:
        os.chdir(TRROSETTA_DIR) 
        comando_predict = [PYTHON_TRROSETTA, "network/predict.py", "-m", "model2019_07", a3m_file, npz_file]
        subprocess.run(comando_predict, check=True)
        print(f"✅ Matriz cruda NPZ guardada en: {npz_file}")
    except Exception as e:
        print(f"\n❌ ERROR: Detalle: {e}")
        sys.exit(1)
    finally:
        os.chdir(BASE_DIR)

    # === PASO 4: EXTRACCIÓN A CSV ===
    print("\n[4/4] Extrayendo matriz de distancias a CSV...")
    try:
        comando_extraer = [PYTHON_TRROSETTA, SCRIPT_EXTRAER, npz_file, csv_file]
        subprocess.run(comando_extraer, check=True)
        print("\n=================================================")
        print(" 🎉 PIPELINE COMPLETADO CON ÉXITO")
        print(f" Matriz final guardada en: {csv_file}")
        print("=================================================")
    except subprocess.CalledProcessError as e:
        print("\n❌ ERROR en el paso 4 (Extracción de CSV). Verifica que el script distance_matrix.trRosetta.py no tenga errores.")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline Automático de trRosetta")
    parser.add_argument("fasta", help="Nombre del archivo FASTA de entrada")
    parser.add_argument("--skip-msa", action="store_true", help="Salta ColabFold y usa el .a3m existente")
    
    args = parser.parse_args()
    ejecutar_pipeline_trRosetta(args.fasta, args.skip_msa)