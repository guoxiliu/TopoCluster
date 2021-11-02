#!/bin/bash

if [ -z $1 ]; then
  printf "Usage: $0  <job_type>\n"
  printf "Please set the job type argument and rerun the program!\n"
  printf "1 for Explicit TopoCluster\n"
  printf "2 for Implicit TopoCluster\n"
  printf "3 for TTK Triangulation\n"
  exit 1
fi

DATASET_FOLDER="../datasets/"

# Job for Explicit TopoCluster
if [ $1 == 1 ]; then
  printf "Testing critical points plugin on Explicit TopoCluster...\n\n"
  RESULT_FILE="results_cp_explicit.txt"
  for dataName in Redsea Engine Cat Sphere Foot Shapes Hole Stent
  do
    printf "\nDataset: ${dataName}\n" &>> "${RESULT_FILE}"
    VTU_FILE="${DATASET_FOLDER}${dataName}_5000.vtu"
    # check if file exists 
    if [ -f $VTU_FILE ]; then 
      scalarFieldCriticalPointsCmd -i "${VTU_FILE}" -F 1 -t 1 &>> "${RESULT_FILE}"
    else 
      printf "Cannot find ${VTU_FILE}! Please double check if the file exists!\n"
    fi
  done

# Job for Implicit TopoCluster
elif [ $1 == 2 ]; then
  printf "Testing critical points plugin on Implicit TopoCluster...\n\n"
  RESULT_FILE="results_cp_implicit.txt"
  for dataName in Redsea Engine Cat Sphere Foot Shapes Hole Stent
  do
    printf "\nDataset: ${dataName}\n" &>> "${RESULT_FILE}"
    VTU_FILE="${DATASET_FOLDER}${dataName}_5000.vtu"
    # check if file exists 
    if [ -f $VTU_FILE ]; then 
      scalarFieldCriticalPointsCmd -i "${VTU_FILE}" -F 1 -t 1 &>> "${RESULT_FILE}"
    else 
      printf "Cannot find ${VTU_FILE}! Please double check if the file exists!\n"
    fi
  done

# Job for TTK Triangulation
elif [ $1 == 3 ]; then
  printf "Testing critical points plugin on TTK Triangulation...\n\n"
  RESULT_FILE="results_cp_ttk.txt"
  DATASET_FOLDER="/home/guoxil/Downloads/"

  for dataName in Redsea Engine Cat Sphere Foot Shapes Hole Stent
  do
    printf "\nDataset: ${dataName}\n" &>> "${RESULT_FILE}"
    VTU_FILE="${DATASET_FOLDER}${dataName}_5000_reindexed.vtu"
    # check if file exists 
    if [ -f $VTU_FILE ]; then 
      scalarFieldCriticalPointsCmd -i "${VTU_FILE}" -t 1 &>> "${RESULT_FILE}"
    else 
      printf "Cannot find ${VTU_FILE}! Please double check if the file exists!\n"
    fi
  done

# Wrong parameter value
else
  printf "Error: Unsupported job type! Please make sure the value is from 1 to 3!\n"
  
fi



