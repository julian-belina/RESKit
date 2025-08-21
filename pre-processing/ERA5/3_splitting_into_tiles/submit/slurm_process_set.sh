#!/bin/bash
#SBATCH --output="logs/slurm-%x-%A-%a.out"
#SBATCH --job-name=ERA5_2015_all_tiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#####SBATCH --array=2010-2018
#####SBATCH --exclude=cn[1-5]

###################################
## set X and Y tiles
TILES_X=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
TILES_Y=(0 1 2 3 4 5 6 7 8 9 10 11 12 13)

declare -a YEARS
declare -a TILES_XX
declare -a TILES_YY

#for yr in $(seq 1980 2019); do

for xi in ${TILES_X[@]}; do
    for yi in ${TILES_Y[@]}; do
        YEARS+=($yr)
        TILES_XX+=($xi)
        TILES_YY+=($yi)
    done
done

#done

#####################################

year=2015 #### ${YEARS[$SLURM_ARRAY_TASK_ID]}
xi=${TILES_XX[$SLURM_ARRAY_TASK_ID]}
yi=${TILES_YY[$SLURM_ARRAY_TASK_ID]}
zoom=4


if [ -z "$year" ]; 
then
    echo "Array id is too high, exiting..."
    exit 0
fi

export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
export MKL_NUM_THREADS=1

source activate cdo

echo "python ./ERA5_processor_tile.py $year $xi $yi $zoom"
python ./ERA5_processor_tile.py $year $xi $yi $zoom

