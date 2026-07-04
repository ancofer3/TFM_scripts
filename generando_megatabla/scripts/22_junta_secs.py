'''
Junta Secs

Fichero: 20_junta_secs.py

Propósito del programa:
Juntamos todas las secuencias que tengan mismo transcript_ID. Los fastas están
ya en orden correcto. Tanto en hebra positiva como negativa, los exones están 
en dirección 5'->3'. 
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def crea_dic(fasta):
    sec_dic = {}
    for record in fasta:
        transcript_id = record.id.split(":")[0]
        if transcript_id not in sec_dic.keys():
            sec_dic[transcript_id] = record.seq
        else:
            sec_dic[transcript_id] += record.seq
    return sec_dic

def saca_records(sec_dic):
    records = []
    for id,sec in sec_dic.items():
        record = SeqRecord(
                        Seq(str(sec)),
                        id=id,
                        description=""
                          )
        records.append(record)
    return records

def saca_fasta(path):
    fasta = SeqIO.parse(path, "fasta")
    sec_dic = crea_dic(fasta)
    records = saca_records(sec_dic)
    SeqIO.write(records,path.replace(".fasta","_fusionados.fasta"),"fasta")

def main():
    for i in ["CDS","3pUTR","5pUTR"]:
        print(i)
        path = f"fastas/secuencias_{i}.fasta"
        saca_fasta(path)

if __name__ == "__main__":
    main()
