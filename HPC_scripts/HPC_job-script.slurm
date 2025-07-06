#!/bin/bash

### All parameters marked with #CHANGE below need to be changed to use this script. Also remove the #CHANGE notations

# Give parameters to queue system: 
#SBATCH --account=<ACCOUNT NAME> #CHANGE
#SBATCH --partition=normal
#SBATCH --job-name=cp4-201028-1plate
#SBATCH --array=1-257%103
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --gres=localscratch:1G
#SBATCH --output=<OUTPUT FILENAME AND LOCATIONS, e.g. .../201028-1_Plate11/slurm-output/slurm-%A_%a.out> #CHANGE

# Turn on all mail notification, and also provide mail address: 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<YOUR EMAIL ADDRESS>


#########################################################################################

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
set -o xtrace   # Get every command printed to stdout

# Set up proper environment:

# Reset the modules to the system default and load Miniconda3/py38_4.8.3
module --quiet purge
module load Miniconda3/py38_4.8.3 #CHANGE depending on version of conda/python installed
module list  # For easier debugging

# Set the ${PS1} (needed in the source of the Anaconda environment) - is this necessary?
export PS1=\$

# Source the conda environment setup
# The variable ${EBROOTMINICONDA3} comes with the module load command
source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh #CHANGE depending on version of conda

# Deactivate any spill-over environment from the login node
conda deactivate &>/dev/null


# Activate the environment by using the full path (not name) 
# to the environment. The full path is listed if you do 
# conda info --envs at the command prompt. 
conda activate <"ENVIRONMENT PATH"> #CHANGE


#########################################################################################

# Move to correct path:
cd <"PLATE PATH, e.g. .../201028-1_Plate11"> #CHANGE

# Find files that need to be copied to localscratch

set +o errexit #Avoid exiting due to error message from cellprofiler

# Batch_data.h5 from CellProfiler must be placed in output-data folder of plate

metadatajsonstring=$(cellprofiler --print-groups ./output-data/Batch_data.h5 2>/dev/null)

set -o errexit #Exit on other errors

JSONINDEX=$(($SLURM_ARRAY_TASK_ID-1))

platedate=$(echo ${metadatajsonstring} | jq -r '.['$JSONINDEX'] | .[0] | .["Metadata_Plate"]')
wellnumber=$(echo ${metadatajsonstring} | jq -r '.['$JSONINDEX'] | .[0] | .["Metadata_Well"]')

#Copy relevant files and make relevant folders

mkdir $LOCALSCRATCH/output-data
cp ./output-data/Batch_data.h5 $LOCALSCRATCH/output-data

#CHANGE - Ensure folders below match your setup
mkdir -p $LOCALSCRATCH/input-images/201028-1_Plate11_rerun1_4645/TimePoint_1 #CHANGE
cp $(find ./input-images/ -name "${platedate}_${wellnumber}_*" -type f) $LOCALSCRATCH/input-images/201028-1_Plate11_rerun1_4645/TimePoint_1 #CHANGE
cp -r ./input-images/sorting-rules $LOCALSCRATCH/input-images #CHANGE
cd $LOCALSCRATCH

cleanup "cp $LOCALSCRATCH/output-data/*.csv <OUTPUT DATA STORAGE PATH>" #CHANGE

# Execute the python program:

# first=$(($SLURM_ARRAY_TASK_ID*4-3))
# last=$(($SLURM_ARRAY_TASK_ID*4))

firstimage=$(echo ${metadatajsonstring} | jq -r '.['$JSONINDEX'] | .[1] | first')
lastimage=$(echo ${metadatajsonstring} | jq -r '.['$JSONINDEX'] | .[1] | last')

python -m cellprofiler -c -r -p ./output-data/Batch_data.h5 -f ${firstimage} -l ${lastimage} 
