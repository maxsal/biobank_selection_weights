#!/bin/bash
#SBATCH -t 12:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 24 ## REQUEST NODES AND CORES (ONLY STAGING NODES CAN REQUEST SINGLE CORES)


mkdir -p TEMP
INPUTFILE="/net/junglebook/home/mmsalva/createUKBphenome/data/UKB-weight-variables.tab"
OUTPUTFILE="/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab"
Rscript --vanilla UKBLoadWeightVariables.R $INPUTFILE
Rscript --vanilla UKBInferAssessmentCentre.R
Rscript --vanilla UKBCreateWeights.R $OUTPUTFILE
