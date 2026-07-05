'''
Sacar propiedades

Fichero: 21_saca_props.py

Propósito del programa:
Recorremos los fastas de 3UTR, 5UTR y CDS y sacamos length y %GC
'''
import glob
import pandas as pd
from Bio import SeqIO


def saca_GC_len(fasta,tipo):
    ids = []
    lengths = []
    gcs = [] 
    for record in fasta:
        transcript_id = record.id
        length = len(record.seq)
        sec = record.seq
        gc = (sec.count("G") + sec.count("C"))/length
        ids.append(transcript_id)
        lengths.append(length)
        gcs.append(gc)
    dic_final = {"transcript_id":ids,
                 f"length_{tipo}":lengths,
                 f"porc_GC_{tipo}":gcs}
    return dic_final

dfs = {}

# Hemos asegurado antes que todos los id que hay en 3pUTR y 5pUTR están en CDS
for i in glob.glob("fastas/*fusionados.fasta"):
    print(i)
    fasta = SeqIO.parse(i,"fasta")
    tipo = i.split("_")[-2]
    transcript_GC = saca_GC_len(fasta,tipo)
    dfs[tipo] = pd.DataFrame(transcript_GC)

# Ahora lo suyo es que creemos el DF entero con los %GC
df_final = dfs["CDS"]
df_final = df_final.merge(dfs["3pUTR"], on="transcript_id", how="left")
df_final = df_final.merge(dfs["5pUTR"], on="transcript_id", how="left")

# Y que lo mergeemos con la megatabla
megatabla = pd.read_csv("megatablas/megatabla_6.tsv",
                        sep="\t")
megatabla_nueva = megatabla.merge(df_final,
                                  left_on="associated_transcript",
                                  right_on="transcript_id",
                                  how="left")
megatabla_nueva = megatabla_nueva.drop(columns=["transcript_id"])
print(megatabla_nueva.isnull().sum())
megatabla_nueva.to_csv("megatablas/megatabla_7.tsv", sep="\t", index=False)
