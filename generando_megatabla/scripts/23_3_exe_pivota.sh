#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=14
#SBATCH -N 1                # Número de nodos
#SBATCH -t 1-00:00:00         # Tiempo máximo (hh:mm:ss)
#SBATCH --qos short         # Cola / QoS
#SBATCH --mem=26G           # Memoria total
#SBATCH -o pivota_%j.out     # Archivo de salida (stdout)
#SBATCH -e pivota_%j.err     # Archivo de errores (stderr)
#SBATCH --job-name=pivota
module load anaconda
conda activate jae_env
python 23_3_pivota.py
