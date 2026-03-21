# Guardar en: ~/Documents/Jazmin/pipeline_esm2.py
import subprocess
import os
import sys

def ejecutar_pipeline_esm2(fasta_filename):
    # 1. Definir rutas
    fasta_path = os.path.join("fasta_files", fasta_filename)
    directorio_resultados = "resultados_esm2"
    
    if not os.path.exists(fasta_path):
        print(f"❌ Error: No se encontró el archivo {fasta_path}")
        sys.exit(1)
        
    # Crear carpeta de resultados si no existe
    if not os.path.exists(directorio_resultados):
        os.makedirs(directorio_resultados)
        
    nombre_base = os.path.splitext(fasta_filename)[0]
    csv_salida = os.path.join(directorio_resultados, f"{nombre_base}_contactos_esm2.csv")
    
    print("=================================================")
    print(f" 🧬 INICIANDO PIPELINE ESM-2 PARA: {nombre_base}")
    print("=================================================")

    # LA MAGIA: Ruta absoluta al Python de tu entorno esmfold
    PYTHON_ESMFOLD = "/home/biocomp/anaconda3/envs/esmfold/bin/python"
    
    # Ruta al script motor que acabamos de guardar en la carpeta esm2/
    SCRIPT_MOTOR = os.path.join("esm2", "generar_matriz_esm2.py")

    try:
        print(f"\nGenerando matriz de contactos con el modelo de lenguaje...")
        print("(Nota: La primera vez descargará el modelo, puede tardar un poco)")
        
        # Armamos el comando: /ruta/a/python esm2/generar_matriz_esm2.py entrada.fasta salida.csv
        comando = [
            PYTHON_ESMFOLD, 
            SCRIPT_MOTOR, 
            fasta_path, 
            csv_salida
        ]
        
        # Ejecutamos!
        subprocess.run(comando, check=True)

        print("\n=================================================")
        print(" 🎉 PIPELINE COMPLETADO CON ÉXITO")
        print(f" Matriz guardada en: {csv_salida}")
        print("=================================================")

    except subprocess.CalledProcessError as e:
        print(f"\n❌ ERROR: El script de ESM-2 falló.")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR INESPERADO: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python pipeline_esm2.py <nombre_archivo.fasta>")
        sys.exit(1)
        
    archivo_fasta_input = sys.argv[1]
    ejecutar_pipeline_esm2(archivo_fasta_input)