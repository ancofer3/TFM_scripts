import pandas as pd
import os
import glob

for fich in glob.glob("kmeros/secuencias*"):

	# Leemos
	#fich = input("Introduce el nombre del fichero: ")
	print(f"Procesando fichero {fich}")
	tipo = os.path.basename(fich).split("_")[1]
	k = os.path.basename(fich).split("_")[-1].split(".")[0]
	df = pd.read_csv(fich, sep="\t")

	# Pivotamos
	pivot = df.pivot(index='SeqID', columns='kmer', values='count').fillna(0)
	pivot['total'] = pivot.sum(axis=1)
	pivot = pivot.div(pivot['total'], axis=0)
	pivot = pivot.drop(["total"], axis=1)
	# Se guarda como como TSV
	pivot = pivot.rename(columns={col:f"kmer_{tipo}_" + col for col in pivot.columns})
	pivot.to_csv(f"pivotadas/{tipo}_{k}mer_pivoted.tsv", sep="\t")
