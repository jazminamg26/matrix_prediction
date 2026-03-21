import torch
import esm
import numpy as np
import pandas as pd  # <-- ¡Novedad añadida!
import sys
import os

def predecir_contactos_esm2(fasta_file, output_csv):
    # 1. Leer el archivo FASTA
    with open(fasta_file, 'r') as f:
        lineas = f.readlines()
        
    secuencia = ""
    for linea in lineas:
        if not linea.startswith(">"):
            secuencia += linea.strip()
            
    L_seq = len(secuencia)
    
    # 2. Cargar el modelo ESM-2 (Modelo de 33 capas)
    print("Cargando modelo ESM-2 (esto puede tomar unos segundos)...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    
    # 3. Mover a GPU si hay una disponible
    if torch.cuda.is_available():
        model = model.cuda()
        
    model.eval()
    
    # 4. Preparar datos
    data = [("proteina", secuencia)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()
        
    # 5. Predecir contactos (Atención)
    print(f"Prediciendo contactos para secuencia de longitud {L_seq}...")
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        
    # 6. Extraer matriz
    contactos = results["contacts"][0].cpu().numpy()
    L_matriz = contactos.shape[0]
    
    # 7. Guardar con formato estricto (headers y sin fila vacía)
    nombres_columnas = [f"c{i+1}" for i in range(L_matriz)]
    df_matriz = pd.DataFrame(contactos, columns=nombres_columnas)
    
    texto_csv = df_matriz.to_csv(index=False, header=True, float_format="%.4f")
    
    # === EL CHEQUEO ESTRICTO para borrar la fila vacía ===
    lineas_csv = texto_csv.split('\n')
    lineas_validas = [linea for linea in lineas_csv if linea.strip() != ""]
    texto_final = '\n'.join(lineas_validas)
    
    with open(output_csv, 'w') as f:
        f.write(texto_final)
        
    print(f"¡Éxito! Matriz ESM-2 ({L_matriz}x{L_matriz}) guardada con encabezados en: {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python generar_matriz_esm2.py <entrada.fasta> <salida.csv>")
        sys.exit(1)
        
    fasta_in = sys.argv[1]
    csv_out = sys.argv[2]
    predecir_contactos_esm2(fasta_in, csv_out)