#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M aarmijo@mit.edu
###############################

###----- RECOMMENDED DEFAULTS -----
picard=/net/bmc-pub7/data0/essigmannlab/dependencies/picard-tools-1.119
GATK=/net/bmc-pub7/data0/essigmannlab/dependencies/GATK
DS_path=/net/bmc-pub7/data0/essigmannlab/jobs/dcs/DCS-3.00
ref=/net/bmc-pub7/data0/essigmannlab/jobs/dcs/ref/EG10_custom.fasta

activate=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin/activate;
ucm=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/envs/ucm/;

endID=''
minMem=3
maxMem=1000
cutOff=0.7
nCutOff=0.3
readLength=133
barcodeLength=8
tagLength=8
spacerLength=5
repFilt=9
readOut=1000000
minDepth=100
minClonality=0
maxClonality=100
softclip_cycles='1-8,120-137'
###----- END -----

###----- NONDEFAULTS -----
###----- END -----

read1=$inp_path/$runID\_$samID\_1\_sequence.fastq
read2=$inp_path/$runID\_$samID\_2\_sequence.fastq

export PATH=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin:$PATH
export PATH=/bin:$PATH

source $activate $ucm;
cd /net/bmc-pub7/data0/essigmannlab/jobs/dcs

date
echo 'Beginning UCM run...'

mkdir $(pwd)/tmp

java -jar -Xmx2g -Djava.io.tmpdir=$(pwd)/tmp $picard/FastqToSam.jar F1=$read1 F2=$read2 O=$run_name.bam SAMPLE_NAME=$run_name TMP_DIR=$(pwd)/tmp

python UnifiedConsensusMaker.py --input $run_name.bam --prefix $run_name --taglen $tagLength --spacerlen $spacerLength --minmem $minMem --maxmem $maxMem --cutoff $cutOff --Ncutoff $nCutOff

gzip -d $run_name'_read1_dcs.fq.gz'
gzip -d $run_name'_read2_dcs.fq.gz'

bwa mem $ref $run_name'_read1_dcs.fq' $run_name'_read2_dcs.fq' > $run_name.dcs.sam

samtools view -Sbu $run_name.dcs.sam | samtools sort -o $run_name.dcs.aln.sort.bam

samtools index $run_name.dcs.aln.sort.bam

samtools view -F 4 -b $run_name.dcs.aln.sort.bam > $run_name.dcs.filt.bam

java -jar $picard/AddOrReplaceReadGroups.jar INPUT=$run_name.dcs.filt.bam OUTPUT=$run_name.dcs.filt.readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

samtools index $run_name.dcs.filt.readgroups.bam

java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $run_name.dcs.filt.readgroups.bam -o $run_name.dcs.filt.readgroups.intervals
java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $run_name.dcs.filt.readgroups.bam -targetIntervals $run_name.dcs.filt.readgroups.intervals -o $run_name.dcs.filt.readgroups.realign.bam
java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T ClipReads -I $run_name.dcs.filt.readgroups.realign.bam -o $run_name.dcs.filt.readgroups.clipped.bam -R $ref --cyclesToTrim $softclip_cycles --clipRepresentation SOFTCLIP_BASES

samtools mpileup -AB -d500000 -f $ref $run_name.dcs.filt.readgroups.clipped.bam > $run_name.dcs.pileup

python $DS_path/CountMuts.py -i $run_name.dcs.pileup -d $minDepth -c $minClonality -C $maxClonality > $run_name.countmuts
python $DS_path/mut-position.py -i $run_name.dcs.pileup -o $run_name.mutpos -d $minDepth -c $minClonality -C $maxClonality
