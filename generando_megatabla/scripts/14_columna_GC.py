'''
------------------------------------------------------------------------
Columna GC

Fichero: 02_columna_GC.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. En este paso empalmamos la tabla con GC creada en el fichero
anterior con la megatabla 

Versión 1
Autores: Andrés Colomer Fernández
Fecha de última modificación: 12/11/2025

------------------------------------------------------------------------
'''
import pandas as pd

megatabla = pd.read_csv("megatablas/megatabla_1.tsv",sep = "\t")
tablaGC = pd.read_csv("tablas_aux/tablaGC.tsv", sep = "\t")

# Por un lado añadimos la columna "exon_junction_density"
megatabla["exon_junction_density"] = (megatabla["ref_exons"]-1)/megatabla["ref_length"]

# por un lado añadimos la tabla GC
megatabla2 = megatabla.merge(
    tablaGC[['RIN', 'TS', 'associated_transcript', 'porc_GC']],
    on=['RIN', 'TS', 'associated_transcript'],
    how='left'
)
print("Valores nulos en df_final:")
print(megatabla2.isnull().sum())
megatabla2.to_csv("megatablas/megatabla_2.tsv",sep="\t",index=False)

