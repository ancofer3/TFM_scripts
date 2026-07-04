#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1 
#SBATCH -t 1-00:00:00         # Tiempo máximo (hh:mm:ss)
#SBATCH --qos short         # Cola / QoS
#SBATCH --mem=100G           # Memoria total
#SBATCH -o logs/random_forest_%j.out     # Archivo de salida (stdout)
#SBATCH -e logs/random_forest_%j.err     # Archivo de errores (stderr)
#SBATCH --job-name=random_forest
# Exportamos las CPUs disponibles para que sklearn las detecte
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load anaconda
conda activate jae_env
python -u 26_random_forest_optuna_v3.py
