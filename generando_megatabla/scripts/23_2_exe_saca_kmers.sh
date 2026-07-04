#!/bin/bash
file=$1
for fich in fastas/*_fusionados.fasta ;do
	for i in {2..7}; do
		scripts/23_1_saca_kmers $fich $i;
		mv ${fich}_K_${i}.txt kmeros/
		done
	done
