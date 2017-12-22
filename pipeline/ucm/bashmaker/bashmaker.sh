#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 1
#$ -M klkim@mit.edu
###############################

run_name=x8901
runID=171017EssA
samID=D17-8901
inp_path=/net/bmc-pub7/data0/essigmannlab/data/171017EssA

# Additional optional parameters
endID=''
ref=''
picard=''
GATK=''
DS_path=''
min=''
max=''
cutOff=''
Ncut=''
rlength=''
blength=''
slength=''
read_type=''
repFilt=''
template=''
softclip_cycles=''

# Generate command
declare -a params=(endID ref picard GATK DS_path min max cutOff Ncut rlength blength slength repFilt template softclip_cycles);

com="python /net/bmc-pub7/data0/essigmannlab/jobs/UCM/bashmaker/UCM_BASH_MAKER.py --run_name $run_name --runIdentifier $runID --sampleIdentifier $samID --inp_path $inp_path"

for p in ${params[@]}; do
  if [ ${!p} ]; then
    com+=" --$p ${!p}"
  fi
done

$com
