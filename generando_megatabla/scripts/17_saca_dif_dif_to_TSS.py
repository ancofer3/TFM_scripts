'''
------------------------------------------------------------------------
Sacar_ diferencia de diferencias al tss

Fichero: 17_saca_dif_dif_to_TSS.py

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
megatabla = pd.read_csv("megatablas/megatabla_4.tsv", sep = '\t')
print(megatabla.columns)
# Nos generamos unas columnas auxiliares que sean 
# El RIN mínimo para una misma comb de transcrito y TS
megatabla["RIN_min"] = megatabla.groupby(["clave","TS"])["RIN"].transform("min")
# El RIN máximo para una misma comb de transcrito y TS
megatabla["RIN_max"] = megatabla.groupby(["clave","TS"])["RIN"].transform("max")
# Vamos a crear subtablas
min = megatabla[megatabla["RIN"] == megatabla["RIN_min"]][["clave",
                                                           "TS",
                                                           "diff_to_gene_TSS"]]
min = min.rename(columns={"diff_to_gene_TSS":"diff_to_gene_TSS_min"})

max = megatabla[megatabla["RIN"] == megatabla["RIN_max"]][["clave",
                                                           "TS",
                                                           "diff_to_gene_TSS"]]
max = max.rename(columns={"diff_to_gene_TSS":"diff_to_gene_TSS_max"})

# Lo unimos a la megatabla de antes
megatabla = megatabla.merge(min, on=['clave', 'TS'], how='left')
megatabla = megatabla.merge(max, on=['clave', 'TS'], how='left')

# Sacamos la diferencia
megatabla['delta(diff_to_gene_TSS)'] = abs(megatabla['diff_to_gene_TSS_min'] - megatabla['diff_to_gene_TSS_max'])

# Quitamos las columnas estas auxiliares
megatabla = megatabla.drop(columns=['RIN_min', 'RIN_max', 'diff_to_gene_TSS_min', 'diff_to_gene_TSS_max'])
print(megatabla.isnull().sum())
megatabla.to_csv(path_or_buf="megatablas/megatabla_5.tsv",sep = "\t", index=False)
