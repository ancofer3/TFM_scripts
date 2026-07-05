'''
------------------------------------------------------------------------
Coverage

Fichero: 18_saca_coverage_prueba.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. Vamos a sacar la cobertura dividida en cuantiles. Empezamos
generando una tabla con las coberturas divididas en cuantiles por 
transcrito para cada TS/RIN. 

La idea es sacar para cada transcrito una lista con el porcentaje de 
reads que caen en cada cuantil.
Para cada input, sacamos cada transcrito unico y de ahi sacamos todos
los diffs_to_gene_TSS. 

Versión 1
Autores: Andrés Colomer Fernández
Fecha de última modificación: 16/04/2026

------------------------------------------------------------------------
'''

import glob as glob 
import pandas as pd
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor

class Classification():
    def __init__(self,path):
        cols_to_use = ['isoform', 'chrom', 'strand', 'length', 'exons', 'structural_category',
                       'associated_gene', 'associated_transcript', 'ref_length', 'ref_exons',
                       'diff_to_gene_TSS', 'diff_to_gene_TTS']
        self.tabla = pd.read_csv(path,sep="\t",usecols=cols_to_use)
    
        # Filtramos los novel
        self.tabla = self.tabla[self.tabla["associated_transcript"] != "novel"]
        fields = os.path.basename(path).split("_")
        self.TS = fields[0]
        self.RIN = fields[2]
        
def divide_and_average(arr, n_fragments=5):
    """Divides a NumPy array into n_fragments and computes the mean of each fragment."""
    fragments = np.array_split(arr, n_fragments)  # Split into nearly equal parts
    means = np.array([np.mean(fragment) for fragment in fragments])  # Compute mean of each
    return means

def sacaCobertura(group):
    # Generamos un vector de 0s del tamaño del transcrito
    res = np.zeros(int(group["ref_length"].iloc[0]))
    # Estas variables son constantes en todas las filas del group
    txStart = 0
    txEnd = group["ref_length"].iloc[0]
    strand = group["strand"].iloc[0]
    # Computamos start y end
    start = (txStart - group['diff_to_gene_TSS']).clip(lower=txStart)
    end = (txEnd + group['diff_to_gene_TTS']).clip(upper=txEnd)
  
    group["start"] = start
    group["end"] = end
    # Corregimos starts y ends 
    group["start"] = group["start"].astype(int)
    group["end"] = group["end"].astype(int)
    
    for _, row in group.iterrows():
        # Para cada fila
        res[int(row["start"]):int(row["end"])+1] += 1
    res = divide_and_average(res,10)
    min_val, max_val = res.min(), res.max()
    if max_val > min_val:
        res = (res-min_val)/(max_val-min_val)
    else:
        res = np.ones_like(res)
    # Para conseguir las columnas con cada decil
    n_bins = len(res)
    col_names = [f"Cov_D{i}" for i in range(1,n_bins+1)]
    return pd.Series(res,index=col_names)


def paralel_saca_cov(path):
    fich = Classification(path)
    # Agrupamos por grupo de asociado
    df_cobertura = fich.tabla.groupby("associated_transcript").apply(sacaCobertura)
    # Añadimos una columna de TS y de RIN
    df_cobertura["TS"] = fich.TS
    df_cobertura["RIN"] = fich.RIN
    return df_cobertura
    

if __name__ == "__main__":
    paths =  glob.glob("inputs/*")
    with ProcessPoolExecutor() as executor:
        resultados = list(executor.map(paralel_saca_cov,paths))
    df_cobertura = pd.concat(resultados)
    df_cobertura.to_csv("tablas_aux/coberturas.tsv",sep= "\t")

