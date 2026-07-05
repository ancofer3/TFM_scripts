'''
Cuenta codones

Fichero: 22_A_cuenta_codones

Propósito del programa:
Sacamos las frecuencias de los codones y la añadimos a una tabla
'''

from Bio import SeqIO
import itertools
import pandas as pd
fasta = SeqIO.parse("fastas/secuencias_CDS_fusionados.fasta", "fasta")
n = 0

def crea_dic_codones():
	codones = itertools.product(['A','C','G','T'], repeat=3)
	codones_count = {}
	for i in codones:
		codones_count["".join(i)]=0
	return codones_count

def cuenta_codones(record):
	codones_count = crea_dic_codones()
	total_codones = 0
	for x in range(0, len(record),3):
		kmero = str(record.seq[x:x+3])
		if kmero in codones_count:
			codones_count[kmero] += 1
			total_codones += 1
	if total_codones > 0:
		for i in codones_count.keys():
			frec_abs = codones_count[i]
			frec_rel = frec_abs/total_codones
			codones_count[i] = frec_rel
	codones_count["transcript_id"] =  record.id.split("::")[0]
	return codones_count

codones = list(crea_dic_codones().keys())
resultados = []
for record in fasta:
	resultados.append(cuenta_codones(record))
df = pd.DataFrame(resultados,columns=["transcript_id"] + codones)
cols = [cod for cod in df.columns if cod != "transcript_id"]
df.rename(columns={cod: f'codon_{cod}' for cod in cols}, inplace=True)
df.to_csv("frec_relativa_codones.tsv", sep="\t", index=False)

