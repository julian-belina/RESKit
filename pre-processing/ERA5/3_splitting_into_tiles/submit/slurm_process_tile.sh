#!/bin/bash
#SBATCH --output="logs/slurm-%x-%A-%a.out"
#SBATCH --job-name=ERA5_Tile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --array=2015
#####SBATCH --exclude=[cn1-5]

year=$SLURM_ARRAY_TASK_ID
xi=$1
yi=$2
zoom=4

export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
export MKL_NUM_THREADS=1

source activate cdo

python ./ERA5_processor_tile_streamlined.py $year $xi $yi $zoom

