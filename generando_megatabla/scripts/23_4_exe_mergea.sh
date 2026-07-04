#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=14
#SBATCH -N 1                # Número de nodos
#SBATCH -t 1-00:00:00         # Tiempo máximo (hh:mm:ss)
#SBATCH --qos short         # Cola / QoS
#SBATCH --mem=1000G           # Memoria total
#SBATCH -o mergea_%j.out     # Archivo de salida (stdout)
#SBATCH -e mergea_%j.err     # Archivo de errores (stderr)
#SBATCH --job-name=mergea
module load anaconda
conda activate jae_env
python -u 23_4_mergea.py
