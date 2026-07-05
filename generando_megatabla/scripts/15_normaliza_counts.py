'''
------------------------------------------------------------------------
Columna GC

Fichero: 15_normaliza_counts.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. En este paso normalizamos los counts por profundidad
a CPM. Sacamos los CPMs a nivel de transcrito y de los subgrupos de 
reads de un mismo transcrito que comparten numero de exones que cubren

Versión 1
Autores: Andrés Colomer Fernández
Fecha de última modificación: 11/02/2026

------------------------------------------------------------------------
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 

tabla = pd.read_csv("megatablas/megatabla_2.tsv", sep= "\t")

# Creamos una subtabla con los valores de RIN-TS y los counts
subtabla = tabla[["RIN","TS","counts","counts_transcript"]]
subtabla["RIN_TS"] = subtabla["RIN"].astype("str") + "_" + subtabla["TS"].astype("str")

subtabla["profundidad"] = subtabla.groupby(by=["RIN","TS"])["counts"].transform("sum")
subtabla["profundidad_counts_transcript"] = subtabla.groupby(by=["RIN","TS"])["counts_transcript"].transform("sum")
subtabla["CPM"] = subtabla["counts"]/(subtabla["profundidad"]/1000000)
subtabla["CPM_transcript"] = subtabla["counts_transcript"]/(subtabla["profundidad_counts_transcript"]/1000000)
tabla["CPM"] = subtabla["CPM"]
tabla["CPM_transcript"] = subtabla["CPM_transcript"]
# Vamos a crear una clave primaria que nos separe todos los transcritos (FSMs y ISMs divididos por exones)
tabla["clave"] = tabla["associated_transcript"] + "_" + tabla["covered_exons"].astype(str) + tabla["structural_category"]
tabla.to_csv("megatablas/megatabla_3.tsv",sep="\t",index=False)
