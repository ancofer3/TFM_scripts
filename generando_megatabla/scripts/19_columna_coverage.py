'''
------------------------------------------------------------------------
Columna Coverage

Fichero: 19_columna_coverage.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. Ya hemos generado una tabla con las coberturas divididas
en cuantiles por transcrito para cada TS/RIN. Ahora toca mergearla a 
la "megatabla"

Versión 1
Autores: Andrés Colomer Fernández
Fecha de última modificación: 26/02/2026

------------------------------------------------------------------------
'''
import pandas as pd
import os
print("FILE:", os.path.abspath(__file__))
megatabla = pd.read_csv("megatablas/megatabla_5.tsv", sep="\t")
coberturas = pd.read_csv("tablas_aux/coberturas.tsv", sep="\t")
nueva = megatabla.merge(coberturas, on=["associated_transcript","TS","RIN"], how="left")
print(nueva.isnull().sum())
nueva.to_csv("megatablas/megatabla_6.tsv", sep="\t", index=False)
