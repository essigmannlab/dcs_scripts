configfile: "config.yaml"

ID = str(list(config["id"].keys())[0])

rule all:
    input:
        expand("post/04/{id}.dcs.aln.sort.bam.d100-c0-0.01.mutpos",id=config["id"])

# Requires earlier processing and run of tag_to_header

###----- PRE-PROCESSING -----

rule index_ref:
    input:
        expand("data/EG10/{reference}.fasta",reference=config["reference"])
    output:
        expand("data/EG10/{reference}.fasta.amb",reference=config["reference"]),
        expand("data/EG10/{reference}.fasta.ann",reference=config["reference"]),
        expand("data/EG10/{reference}.fasta.bwt",reference=config["reference"]),
        expand("data/EG10/{reference}.fasta.pac",reference=config["reference"]),
        expand("data/EG10/{reference}.fasta.sa",reference=config["reference"])
    shell:
        "bwa index {input}"

rule align_to_ref:
    input:
        ref = expand("data/EG10/{reference}.fasta",reference=config["reference"]),
        read = expand("pipeline/01/{sample}.fq.smi",sample=config["samples"])
    output:
        expand("pipeline/02/{sample}.aln",sample=config["samples"])
    run:
        for x in config["samples"]:
            shell("bwa aln {input.ref} pipeline/01/"+x+".fq.smi > pipeline/02/" +x+".aln")

rule make_sam:
    input:
        ref = expand("data/EG10/{reference}.fasta",reference=config["reference"]),
        smi = expand("pipeline/01/{sample}.fq.smi",sample=config["samples"]),
        aln = expand("pipeline/02/{sample}.aln",sample=config["samples"])
    output:
        expand("pipeline/03/{id}.pe.sam",id=config["id"])
    shell:
        "bwa sampe -s {input.ref} {input.aln} {input.smi} > {output}"

rule bam_sort01:
    input:
        expand("pipeline/03/{id}.pe.sam",id=config["id"])
    output:
        expand("pipeline/04/{id}.pe.sort.bam",id=config["id"])
    run:
        for id in config["id"]:
            shell("samtools view -Sbu {input} | samtools sort - pipeline/04/"+str(id)+".pe.sort")

rule consensus_maker:
    input:
        id = expand("pipeline/04/{id}.pe.sort.bam",id=config["id"])
    output:
        tag = expand("pipeline/05/{id}.pe.tagcounts",id=config["id"]),
        bam = expand("pipeline/05/{id}.sscs.bam",id=config["id"]),
        nm = expand("pipeline/05/{id}.sscs_NM.bam",id=config["id"]),
        lcc = expand("pipeline/05/{id}.sscs_LCC.bam",id=config["id"]),
        stats = expand("pipeline/05/{id}.pe.tagstats",id=config["id"])
    shell: # ideally, other inputs would be coded into config
        "python DCS-2.00/ConsensusMaker.py --infile {input.id} "
            "--tagfile {output.tag} "
            "--outfile {output.bam} "
            "--min 3 "
            "--max 1000 "
            "--cutoff 0.7 "
            "--Ncutoff 0.3 "
            "--readlength 133 "
#            "--readlength 84 "
            "--read_type dpm "
            "--filt osn"

rule bam_sort02:
    input:
        expand("pipeline/05/{id}.sscs.bam",id=config["id"])
    output:
        expand("pipeline/06/{id}.sscs.sort.bam",id=config["id"])
    run:
        for id in config["id"]:
            shell("samtools view -bu {input} | samtools sort - pipeline/06/"+str(id)+".sscs.sort")

rule make_dcs:
    input:
        expand("pipeline/06/{id}.sscs.sort.bam",id=config["id"])
    output:
        bam = expand("pipeline/07/{id}.dcs.bam",id=config["id"]),
        r1 = expand("pipeline/07/{id}.dcs.r1.fq",id=config["id"]),
        r2 = expand("pipeline/07/{id}.dcs.r2.fq",id=config["id"])
    shell: # again, ideally inputs would be coded into config
        "python DCS-2.00/DuplexMaker.py --infile {input} --outfile {output.bam} --Ncutoff 0.3 --readlength 133"

rule align_dcs01:
    input:
        ref = expand("data/EG10/{reference}.fasta",reference=config["reference"]),
        r1 = expand("pipeline/07/{id}.dcs.r1.fq",id=config["id"]),
        r2 = expand("pipeline/07/{id}.dcs.r2.fq",id=config["id"])
    output:
        r1 = expand("pipeline/08/{id}.dcs.r1.aln",id=config["id"]),
        r2 = expand("pipeline/08/{id}.dcs.r2.aln",id=config["id"])
    run:
        shell("bwa aln {input.ref} {input.r1} > {output.r1}"),
        shell("bwa aln {input.ref} {input.r2} > {output.r2}")

rule align_dcs02:
    input:
        ref = expand("data/EG10/{reference}.fasta",reference=config["reference"]),
        fr1 = expand("pipeline/07/{id}.dcs.r1.fq",id=config["id"]),
        fr2 = expand("pipeline/07/{id}.dcs.r2.fq",id=config["id"]),
        ar1 = expand("pipeline/08/{id}.dcs.r1.aln",id=config["id"]),
        ar2 = expand("pipeline/08/{id}.dcs.r2.aln",id=config["id"])
    output:
        expand("pipeline/08/{id}.dcs.sam",id=config["id"])
    shell:
        "bwa sampe -s {input.ref} {input.ar1} {input.ar2} {input.fr1} {input.fr2} > {output}"

rule sort_aligned_dcs:
    input:
        expand("pipeline/08/{id}.dcs.sam",id=config["id"])
    output:
        expand("pipeline/09/{id}.dcs.aln.sort.bam",id=config["id"])
    run:
        shell("samtools view -Sbu {input} | samtools sort - pipeline/09/"+ID+".dcs.aln.sort")

rule index_sorted_dcs:
    input:
        expand("pipeline/09/{id}.dcs.aln.sort.bam",id=config["id"])
    output:
        expand("pipeline/09/{id}.dcs.aln.sort.bam.bai",id=config["id"])
    shell:
        "samtools index {input}"

###----- POST-PROCESSING -----

rule filter_01:
    input:
        expand("pipeline/09/{id}.dcs.aln.sort.bam",id=config["id"])
    output:
        expand("post/01/{id}.dcs.filt.bam",id=config["id"])
    shell:
        "samtools view -F 4 -b {input} > {output}"

rule prep_fasta:
    input:
        ref = expand("data/{reference}/{reference}.fasta",reference=config["reference"])
    output:
        d = expand("data/{reference}/{reference}.dict",reference=config["reference"]),
        fai = expand("data/{reference}/{reference}.fasta.fai",reference=config["reference"])
    run:
        shell("java -jar ~/dependencies/picard-tools-1.119/CreateSequenceDictionary.jar R= {input.ref} O= {output.d}"),
        shell("samtools faidx {input.ref}")

rule clip_file:
    input:
        expand("post/01/{id}.dcs.filt.bam",id=config["id"])
    output:
        expand("post/02/{id}.dcs.filt.readgroups.bam",id=config["id"])
    run:
        shell("java -jar ~/dependencies/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT={input} OUTPUT={output} RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default"),
        shell("samtools index {output}")

rule filter_02:
    input:
        b = expand("post/02/{id}.dcs.filt.readgroups.bam",id=config["id"]),
        ref = expand("data/{reference}/{reference}.fasta",reference=config["reference"])
    output:
        i = expand("post/03/{id}.dcs.filt.readgroups.intervals",id=config["id"]),
        r = expand("post/03/{id}.dcs.filt.readgroups.realign.bam",id=config["id"]),
        c = expand("post/03/{id}.dcs.filt.readgroups.clipped.bam",id=config["id"])
    run:
        shell("java -Xmx2g -jar ~/dependencies/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {input.ref} -I {input.b} -o {output.i}"),
        shell("java -Xmx2g -jar ~/dependencies/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R {input.ref} -I {input.b} -targetIntervals {output.i} -o {output.r}"),
        shell("java -Xmx2g -jar ~/dependencies/GATK/GenomeAnalysisTK.jar -T ClipReads -I {output.r} -o {output.c} -R {input.ref} --cyclesToTrim '1-4,117-137' --clipRepresentation SOFTCLIP_BASES")

rule pileup:
    input:
        b = expand("post/03/{id}.dcs.filt.readgroups.clipped.bam",id=config["id"]),
        ref = expand("data/{reference}/{reference}.fasta",reference=config["reference"])
    output:
        p = expand("post/04/{id}.dcs.pileup",id=config["id"]),
    run:
        shell("samtools mpileup -B -A -d 500000 -f {input.ref} {input.b} > {output.p}"),

rule countmuts:
    input:
        b = expand("post/03/{id}.dcs.filt.readgroups.clipped.bam",id=config["id"]),
        p = expand("post/04/{id}.dcs.pileup",id=config["id"]),
        ref = expand("data/{reference}/{reference}.fasta",reference=config["reference"])
    output:
        c = expand("post/04/{id}.dcs.aln.sort.bam.d100-c0-0.01.unique.countmuts",id=config["id"]),
        m = expand("post/04/{id}.dcs.aln.sort.bam.d100-c0-0.01.mutpos",id=config["id"])
    run:
        shell("python DCS-2.00/CountMuts.py -i {input.p} -d 100 -c 0 -C 0.01 > {output.c}"),
        shell("python DCS-2.00/mut-position.py -i {input.p} -o {output.m} -d 100 -c 0 -C 0.01")
