#!/bin/bash

#SBATCH --job-name=dettetore6110
#SBATCH --array=1-26247%40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=%j.%a.o
#SBATCH --error=%j.%a.e
#SBATCH --qos=1day
#SBATCH --mail-user=crstp.strt@gmail.com
#SBATCH --mail-type=END,FAIL

## Adapt the job array to the number of files you have!
## Running this script took ...

GNUMBERS=/scicore/home/gagneux/stritt0001/TB/projects/detettoreTB/NOTEBOOKS/D_SR_data/GENOMES_selection.published_lineage-assigned.max5000.g_number.txt

V2=/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/v2
detettore='/scicore/home/gagneux/stritt0001/programs/detettore6110/detettore6110.py'

source ~/miniconda3/etc/profile.d/conda.sh
conda activate detettore6110

GNR=$(head -n $SLURM_ARRAY_TASK_ID $GNUMBERS | tail -n 1 | cut -f1)
CRAM=${V2}/${GNR:0:3}/${GNR:3:2}/${GNR:5:2}/${GNR}.cram

${detettore} ${CRAM}

mv $SLURM_JOBID.$SLURM_ARRAY_TASK_ID.e ${GNR}.e
mv $SLURM_JOBID.$SLURM_ARRAY_TASK_ID.o ${GNR}.o
