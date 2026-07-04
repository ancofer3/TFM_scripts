'''
Filtra GTF

Fichero: 20_filtra_gtf.py

Propósito del programa:
Filtramos el gtf de anotación para quedarnos con solo las filas que 
corresponden a transcritos que estén en nuestra tabla.
'''
import pandas as pd
import csv

fichero = "/storage/gge/genomes/human/Homo_sapiens.GRCh38.99.gtf"
columnas = ["chr","source","feature","start","end","score","strand","frame","attributes"]
gtf = pd.read_csv(fichero, sep = "\t",comment="#",header=None,names=columnas)
megatabla = pd.read_csv("megatablas/megatabla_6.tsv", sep="\t")
transcritos_unicos = megatabla["associated_transcript"].unique().tolist()
gtf["transcript_id"] = gtf["attributes"].str.extract(r'transcript_id "([^"]+)"')
gtf = gtf[gtf["transcript_id"].isin(transcritos_unicos)]
gtf = gtf.drop(columns="transcript_id")
gtf.to_csv("tablas_aux/anot_filtrada.gtf", sep="\t", quoting=csv.QUOTE_NONE, escapechar="\\", index=False)

