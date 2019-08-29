

###################
# FLAG definition #
###################

# Request the Bourne Shell
#$ -S /bin/bash

# Change to the current working directory before starting the job
#$ -cwd

# Change the job name to "hello_world"

#$ -N TTbar_13TeV_TuneCUETP8M1_cfi_GEN_SIM_10_chunk0 

# Resource request. We request 1MB of memory, and 60 seconds of wall
# clock time, that more than is enough for the test.
# -l mem_free=2G
# -l h_rt=4:0:00

# We are using the "l.tests" project for the examples
#$ -P l.gaes


#################
# Actual script #
#################

pushd "ABSOLUTE PATH WHERE THE CMSSW/SRC DIR RELEASE IS LOCATED"
eval `scramv1 runtime -sh`
pushd
python "ABSOLUTE PATH TO plotDistributions.py running file" 
