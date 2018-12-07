#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 1
#$ -M aarmijo@mit.edu
###############################

#----- Parameters to change -----
run_name=x8897
runID=171017EssA
samID=D17-8897
inp_path=/net/bmc-pub7/data0/essigmannlab/data/$runID

# If FASTQ files include flowcell ID in the file name, change
# otherwise leave endID=''
endID=2833L

#----- Additional optional parameters -----
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

cd /net/bmc-pub7/data0/essigmannlab/jobs/dcs

# Generate command
declare -a params=(endID ref picard GATK DS_path min max cutOff Ncut rlength blength slength repFilt template softclip_cycles);

com="python /net/bmc-pub7/data0/essigmannlab/jobs/dcs/bashmaker/UCM_BASH_MAKER.py --run_name $run_name --runIdentifier $runID --sampleIdentifier $samID --inp_path $inp_path"

for p in ${params[@]}; do
  if [ ${!p} ]; then
    com+=" --$p ${!p}"
  fi
done

$com
