#!/bin/bash

# Filtramos el gtf a solo las regiones CDS, 3pUTR y 5pUTR
gtf="tablas_aux/anot_filtrada.gtf"
# Para la CDS
echo "Sacando la anot de CDS"
awk '$3 == "CDS" { 
    if (match($0, /transcript_id "([^"]+)"/, a)) {
        print $1 "\t" ($4-1) "\t" $5 "\t" a[1] "\t" "." "\t" $7 
    }
}' $gtf > tablas_aux/anotacion_CDS.bed

# Para la 3' UTR
echo "Sacando la anot de 3pUTR"
awk '$3 == "three_prime_utr" { 
    if (match($0, /transcript_id "([^"]+)"/, a)) {
        print $1 "\t" ($4-1) "\t" $5 "\t" a[1] "\t" "." "\t" $7 
    }
}' $gtf > tablas_aux/anotacion_3pUTR.bed

# Para la 5pUTR
echo "Sacando la anot de 5pUTR"
awk '$3 == "five_prime_utr" { 
    if (match($0, /transcript_id "([^"]+)"/, a)) {
        print $1 "\t" ($4-1) "\t" $5 "\t" a[1] "\t" "." "\t" $7 
    }
}' $gtf > tablas_aux/anotacion_5pUTR.bed
# Ahora sacamos las secuencias asociadas
assembly="/storage/gge/genomes/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
module load bedtools
echo "Sacando el fasta de CDS"
bedtools getfasta -fi $assembly -bed tablas_aux/anotacion_CDS.bed -name -s -fo fastas/secuencias_CDS.fasta
echo "Sacando el fasta de 3pUTR"
bedtools getfasta -fi $assembly -bed tablas_aux/anotacion_3pUTR.bed -name -s -fo fastas/secuencias_3pUTR.fasta
echo "Sacando el fasta de 5pUTR"
bedtools getfasta -fi $assembly -bed tablas_aux/anotacion_5pUTR.bed -name -s -fo fastas/secuencias_5pUTR.fasta
echo "FIN"
