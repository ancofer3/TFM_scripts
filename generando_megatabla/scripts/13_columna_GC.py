'''
------------------------------------------------------------------------
Columna GC : 
Fichero: 13_columna_GC.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. En este paso sacamos la columna de porcentajes de GC. Podemos
hacer algo como esto:

creamos una tabla auxiliar 
for muestra in muestras:
    PROCESAMOS EL BAM en un dict id_secuencia
    proceso el input o la megatabla filtrando por ts y rin
    for transcrito unico in tabla:
        secs = []
        for id asociado al transcrito:
            Sacamos la secuencia del sam
            secs.append(secuencia)
        %GC = cuentaGC(secs)
        Añadimos %GC a la tabla con TS, RIN y transcrito

Versión 2
Autores: Andrés Colomer Fernández
Fecha de última modificación: 11/04/2026

------------------------------------------------------------------------
'''
import pysam
import pandas as pd
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def procesa_Bam(path):
    bam = pysam.AlignmentFile(path,"rb")
    secuencias = {}
    for read in bam:
        secuencias[read.query_name] = read.query_sequence
    return secuencias

def sacaGC(secuencias):
    porcentajes = []
    for sec in secuencias:
        porc = 100 * ((sec.count("G") + sec.count("C")) / len(sec))
        porcentajes.append(porc)           
    porc_GC = np.mean(porcentajes)
    return porc_GC

def procesaFichero(fichero):
    filas = []
    # Nos quedamos con el rin y el ts del fichero
    campos = fichero.split("/")[1].split("_")
    print(campos)
    ts = campos[0]
    rin = campos[2]
    # La ruta al bam asociado con las secuencias
    ruta_bam = f"bams_dRNA/{ts}_RIN_{rin}_pass_primary_aln_sorted.bam"
    # Sacamos diccionario id_secuencia
    id_secuencias = procesa_Bam(ruta_bam)
    # Leemos la tabla input como un df
    tabla = pd.read_csv(fichero, sep = "\t", usecols=["isoform", "associated_transcript"])
    # Para cada transcrito único de la tabla
    grupos = tabla.groupby("associated_transcript")["isoform"]
    for transcript, ids_series in grupos:
        # Tomamos todas los ids asociados a ese transcrito
        ids = ids_series.to_list()
        # Y sacamos con el diccionario las secuencia de interés 
        secs = [id_secuencias[id] for id in ids]
        # Sacamos el porcentaje medio de GC
        porc_GC = sacaGC(secs)
        # Y nos lo guardamos todo en una lista para luego meter en una tabla auxiliar
        fila = [ts,rin,transcript,porc_GC]
        filas.append(fila)
    return filas


if __name__ == "__main__":

    ficheros = glob.glob("inputs/*")

    lista_filas = []
    with ProcessPoolExecutor(max_workers=15) as executor:
        for filas in executor.map(procesaFichero, ficheros):
            lista_filas.extend(filas)

    df_final = pd.DataFrame(lista_filas, columns=["TS", "RIN", "associated_transcript", "porc_GC"])
    df_final.to_csv("tablas_aux/tablaGC.tsv", sep="\t", index=False)


