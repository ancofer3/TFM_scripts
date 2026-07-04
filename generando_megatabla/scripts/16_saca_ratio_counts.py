'''
------------------------------------------------------------------------
Sacar ratio_counts

Fichero: 16_saca_ratio_counts.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. En este paso sacamos un ratio de counts en RIN min vs RIN max
para cada transcrito

Versión 1
Autores: Andrés Colomer Fernández
Fecha de última modificación: 11/02/2026

------------------------------------------------------------------------
'''
import pandas as pd
import glob

def SacaRINyTS(path):
    fields = path.split("\\")[-1].split("_")
    TS = fields[0]
    RIN = fields[2]
    return RIN, TS

# Un diccionario que contenga el RIN max y min para cada TS
TS_RINrange = {}
for i in glob.glob("inputs/*"):
    RIN, TS = SacaRINyTS(i)
    RIN = float(RIN)
    if TS not in TS_RINrange:
        TS_RINrange[TS] = [RIN,RIN]
    else:
        min_actual = TS_RINrange[TS][0]
        max_actual = TS_RINrange[TS][1]
        TS_RINrange[TS] = [min(min_actual, RIN), max(max_actual, RIN)]
    rango = TS_RINrange[TS]
    
print(TS_RINrange)
megatabla = pd.read_csv("megatablas/megatabla_3.tsv", sep = '\t')
print(megatabla.columns)
# Nos generamos unas columnas auxiliares que sean 
# El RIN mínimo para una misma comb de transcrito y TS
megatabla["RIN_min"] = megatabla.groupby(["clave","TS"])["RIN"].transform("min")
# El RIN máximo para una misma comb de transcrito,TS y  covered exons
megatabla["RIN_max"] = megatabla.groupby(["clave","TS"])["RIN"].transform("max")
# Vamos a crear subtablas
min = megatabla[megatabla["RIN"] == megatabla["RIN_min"]][["clave",
                                                           "TS",
                                                           "CPM"]]
min = min.rename(columns={"CPM":"counts_min"})

max = megatabla[megatabla["RIN"] == megatabla["RIN_max"]][["clave",
                                                           "TS",
                                                           "CPM"]]
max = max.rename(columns={"CPM":"counts_max"})

# Lo unimos a la megagtabla de antes
megatabla = megatabla.merge(min, on=['clave', 'TS'], how='left')
megatabla = megatabla.merge(max, on=['clave', 'TS'], how='left')

# Sacamos el ratio
megatabla['ratio_counts'] = megatabla['counts_min'] / megatabla['counts_max']

# Quitamos las columnas estas auxiliares
megatabla = megatabla.drop(columns=['RIN_min', 'RIN_max', 'counts_min', 'counts_max'])
print(megatabla.isnull().sum())
megatabla.to_csv(path_or_buf="megatablas/megatabla_4.tsv",sep = "\t", index=False)
