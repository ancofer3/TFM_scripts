#!/bin/bash
ENV_LOAD="module load anaconda && conda activate jae_env"
# ==============================================================================
# PASO 1: Filtrado de classifications)
# ==============================================================================
JOB_1=$(sbatch --parsable \
    --job-name="filtrado_classifications" \
    --qos=short \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso1_filter_%j.out" \
    --error="logs/paso1_filter_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/12_prepara_transcritos_megatabla.py")

echo "Paso 1 enviado con ID: $JOB_1 (1 CPUs, 32GB RAM)"

JOB_2=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_1 \
    --job-name="crea_tabla_GC" \
    --cpus-per-task=15 \
    --mem=64G \
    --time=01:00:00 \
    --output="logs/paso2_GC_%j.out" \
    --error="logs/paso2_GC_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/13_columna_GC.py")

echo "Paso 2 enviado con ID: $JOB_2 (15 CPUs, 64GB RAM)"

JOB_3=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_2 \
    --job-name="crea_col_GC" \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso3_GC_%j.out" \
    --error="logs/paso3_GC_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/14_columna_GC.py")

echo "Paso 3 enviado con ID: $JOB_3 (1 CPUs, 32GB RAM)"

JOB_4=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_3 \
    --job-name="normaliza_counts"  \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso4_CPM_%j.out" \
    --error="logs/paso4_CPM_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/15_normaliza_counts.py")

echo "Paso 4 enviado con ID: $JOB_4 (1 CPUs, 32GB RAM)"

JOB_5=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_4 \
    --job-name="ratio_counts"  \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso5_ratio_%j.out" \
    --error="logs/paso5_ratio_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/16_saca_ratio_counts.py")

echo "Paso 5 enviado con ID: $JOB_5 (1 CPUs, 32GB RAM)"

JOB_6=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_5 \
    --job-name="deltadiff"  \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso6_delta_%j.out" \
    --error="logs/paso6_delta_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/17_saca_dif_dif_to_TSS.py")

echo "Paso 6 enviado con ID: $JOB_6 (1 CPUs, 32GB RAM)"

JOB_7=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_6 \
    --job-name="tabla_cov"  \
    --cpus-per-task=15 \
    --mem=64G \
    --time=01:00:00 \
    --output="logs/paso7_tablaCov_%j.out" \
    --error="logs/paso7_tablaCov_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/18_saca_coverage_prueba.py")

echo "Paso 7 enviado con ID: $JOB_7 (15 CPUs, 64GB RAM)"

JOB_8=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_7 \
    --job-name="columnas_cov"  \
    --cpus-per-task=1 \
    --mem=32G \
    --time=01:00:00 \
    --output="logs/paso8_colCov_%j.out" \
    --error="logs/paso8_colCov_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/19_columna_coverage.py")

echo "Paso 8 enviado con ID: $JOB_8 (1 CPUs, 32GB RAM)"

JOB_9=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_8 \
    --job-name="filtra_gtf" \
    --cpus-per-task=1 \
    --mem=64G \
    --time=01:00:00 \
    --output="logs/paso9_filtraGtf_%j.out" \
    --error="logs/paso9_filtraGtf_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/20_filtra_gtf.py")

echo "Paso 9 enviado con ID: $JOB_9 (1 CPUs, 32GB RAM)"

JOB_10=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_9 \
    --job-name="saca_secs" \
    --cpus-per-task=15 \
    --mem=26G \
    --time=01:00:00 \
    --output="logs/paso10_sacaSecs_%j.out" \
    --error="logs/paso10_sacaSecs_%j.err" \
    --wrap="$ENV_LOAD && bash scripts/21_saca_secs.sh")

echo "Paso 10 enviado con ID: $JOB_10 (15 CPUs, 26GB RAM)"

JOB_11=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_10 \
    --job-name="junta_secs" \
    --cpus-per-task=1 \
    --mem=64G \
    --time=01:00:00 \
    --output="logs/paso11_juntaSecs_%j.out" \
    --error="logs/paso11_juntaSecs_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/22_junta_secs.py")

echo "Paso 11 enviado con ID: $JOB_11 (15 CPUs, 64GB RAM)"

JOB_12=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_11 \
    --job-name="saca_props" \
    --cpus-per-task=1 \
    --mem=64G \
    --time=01:00:00 \
    --output="logs/paso12_sacaProps_%j.out" \
    --error="logs/paso12_sacaProps_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/22_Z_saca_props.py")

echo "Paso 12 enviado con ID: $JOB_12 (1 CPUs, 32GB RAM)"

JOB_13=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_12 \
    --job-name="saca_codones" \
    --cpus-per-task=14 \
    --mem=26G \
    --time=01:00:00 \
    --output="logs/paso13_cuentaCodones_%j.out" \
    --error="logs/paso13_cuentaCodones_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/22_A_cuenta_codones.py")

echo "Paso 13 enviado con ID: $JOB_13 (14 CPUs, 26GB RAM)"

JOB_14=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_13 \
    --job-name="columnas_codones" \
    --cpus-per-task=14 \
    --mem=26G \
    --time=01:00:00 \
    --output="logs/paso14_columnasCodones_%j.out" \
    --error="logs/paso14_columnasCodones_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/22_B_columnas_codones.py")

echo "Paso 14 enviado con ID: $JOB_14 (14 CPUs, 26GB RAM)"

JOB_15=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_14 \
    --job-name="saca_kmers" \
    --cpus-per-task=1 \
    --mem=26G \
    --time=01:00:00 \
    --output="logs/paso15_sacakmers_%j.out" \
    --error="logs/paso15_sacakmers_%j.err" \
    --wrap="$ENV_LOAD && bash scripts/23_2_exe_saca_kmers.sh")
echo "Paso 15 enviado con ID: $JOB_15 (1 CPUs, 26GB RAM)"

JOB_16=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_15 \
    --job-name="pivota" \
    --cpus-per-task=14 \
    --mem=26G \
    --time=01:00:00 \
    --output="logs/paso16_pivota_%j.out" \
    --error="logs/paso16_pivota_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/23_3_pivota.py")

echo "Paso 16 enviado con ID: $JOB_16 (14 CPUs, 26GB RAM)"

JOB_17=$(sbatch --parsable \
    --qos=short \
    --dependency=afterok:$JOB_16 \
    --job-name="mergea" \
    --cpus-per-task=14 \
    --mem=1000G \
    --time=01:00:00 \
    --output="logs/paso17_mergea_%j.out" \
    --error="logs/paso17_mergea_%j.err" \
    --wrap="$ENV_LOAD && python -u scripts/23_4_mergea.py")

echo "Paso 17 enviado con ID: $JOB_16 (14 CPUs, 1000GB RAM)"
