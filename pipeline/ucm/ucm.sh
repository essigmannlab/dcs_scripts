#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 4
#$ -M klkim@mit.edu
########################

# Parameters by run
run_name=x8901
read1=/net/bmc-pub7/data0/essigmannlab/data/171017EssA/171017EssA_D17-8901_1_sequence.fastq
read2=/net/bmc-pub7/data0/essigmannlab/data/171017EssA/171017EssA_D17-8901_2_sequence.fastq

# Generally constant
ref=/net/bmc-pub7/data0/essigmannlab/jobs/UCM/ref/EG10_custom.fasta
picard=/net/bmc-pub7/data0/essigmannlab/dependencies/picard-tools-1.119
GATK=/net/bmc-pub7/data0/essigmannlab/dependencies/GATK
DSpath=/net/bmc-pub7/data0/essigmannlab/jobs/UCM/DCS-3.00

# Environment setup
activate=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin/activate;
ucm=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/envs/ucm/;

export PATH=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin:$PATH
export PATH=/bin:$PATH

source $activate $ucm;

# Begin run
mkdir /net/bmc-pub7/data0/essigmannlab/jobs/UCM/tmp

java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp $picard/FastqToSam.jar F1=$read1 F2=$read2 O=$run_name.bam SAMPLE_NAME=$run_name TMP_DIR=`pwd`/tmp

python UnifiedConsensusMaker.py --input $run_name.bam --prefix $run_name --taglen 8 --spacerlen 5 --minmem 3 --maxmem 1000 --cutoff 0.7 --Ncutoff 0.3

gzip -d $run_name'_read1_dcs.fq.gz'
gzip -d $run_name'_read2_dcs.fq.gz'

bwa mem $ref $run_name'_read1_dcs.fq' $run_name'_read2_dcs.fq' > $run_name.dcs.sam

#bwa aln $ref $run_name'_read1_dcs.fq' > $run_name'_read1_dcs.aln'
#bwa aln $ref $run_name'_read2_dcs.fq' > $run_name'_read2_dcs.aln'
#bwa sampe -s $ref $run_name'_read1_dcs.aln' $run_name'_read2_dcs.aln' $run_name'_read1_dcs.fq' $run_name'_read2_dcs.fq' > $run_name.dcs.sam

samtools view -Sbu $run_name.dcs.sam | samtools sort -o $run_name.dcs.aln.sort.bam

samtools index $run_name.dcs.aln.sort.bam

samtools view -F 4 -b $run_name.dcs.aln.sort.bam > $run_name.dcs.filt.bam

java -jar $picard/AddOrReplaceReadGroups.jar INPUT=$run_name.dcs.filt.bam OUTPUT=$run_name.dcs.filt.readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

samtools index $run_name.dcs.filt.readgroups.bam

java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $run_name.dcs.filt.readgroups.bam -o $run_name.dcs.filt.readgroups.intervals
java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $run_name.dcs.filt.readgroups.bam -targetIntervals $run_name.dcs.filt.readgroups.intervals -o $run_name.dcs.filt.readgroups.realign.bam
java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T ClipReads -I $run_name.dcs.filt.readgroups.realign.bam -o $run_name.dcs.filt.readgroups.clipped.bam -R $ref --cyclesToTrim '1-8,120-137' --clipRepresentation SOFTCLIP_BASES

samtools mpileup -AB -d500000 -f $ref $run_name.dcs.filt.readgroups.clipped.bam > $run_name.dcs.pileup

python $DSpath/CountMuts.py -i $run_name.dcs.pileup -d 100 -c 0 -C 100 > $run_name.countmuts
python $DSpath/mut-position.py -i $run_name.dcs.pileup -o $run_name.mutpos -d 100 -c 0 -C 100
