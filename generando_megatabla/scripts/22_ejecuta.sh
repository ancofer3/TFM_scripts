#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=14
#SBATCH -N 1                # Número de nodos
#SBATCH -t 1-00:00:00         # Tiempo máximo (hh:mm:ss)
#SBATCH --qos short         # Cola / QoS
#SBATCH --mem=26G           # Memoria total
#SBATCH -o columnas_codones_%j.out     # Archivo de salida (stdout)
#SBATCH -e columnas_codones_%j.err     # Archivo de errores (stderr)
#SBATCH --job-name=columnas_codones

module load anaconda
conda activate jae_env
python 22_cuenta_codones.py
python 22_columnas_codones.py
