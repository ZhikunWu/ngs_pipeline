rule fastp:
    """
    Do not trim partial sequence with low quality.
    Discard whole read if the low quality sequence percentage is lower than the threshold.
    """
    input:
        # R1 = IN_PATH + '/raw/{sample}.R1.fastq.gz',
        # R2 = IN_PATH + '/raw/{sample}.R2.fastq.gz',
        R1 = IN_PATH + '/raw/{sample}_1.fastq.gz',
        R2 = IN_PATH + '/raw/{sample}_2.fastq.gz',
    output:
        R1 = IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.clean.paired.R2.fq.gz',
    threads:
        THREADS
    # params:
    #     length_required = config["length_required"],
    #     qualified_quality_phred = config["qualified_quality_phred"],
    #     unqualified_percent_limit = config["unqualified_percent_limit"],
    log:
        IN_PATH + "/log/trim/{sample}_fastp.log", 
    run:
        ### --disable_adapter_trimming --length_required {params.length_required}  --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit}
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 2  >{log} 2>&1")


rule fastqc:
    input:
        raw_R1 = IN_PATH + '/raw/{sample}.R1.fastq.gz',
        raw_R2 = IN_PATH + '/raw/{sample}.R2.fastq.gz',
        clean_R1 = IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz',
        clean_R2 = IN_PATH + '/clean/{sample}.clean.paired.R2.fq.gz',
    output:
        raw_R1 = IN_PATH + '/FastQC/raw/{sample}/{sample}.R1_fastqc.zip', 
        raw_R2 = IN_PATH + '/FastQC/raw/{sample}/{sample}.R2_fastqc.zip',
        clean_R1 = IN_PATH + '/FastQC/clean/{sample}/{sample}.clean.paired.R1_fastqc.zip',  
        clean_R2 = IN_PATH + '/FastQC/clean/{sample}/{sample}.clean.paired.R2_fastqc.zip',
    threads:
        THREADS
    log:
        raw = IN_PATH + "/log/FastQC/{sample}_raw.log",
        clean = IN_PATH + "/log/FastQC/{sample}_clean.log",
    params:
        FASTQC = BIN_DIR + '/FastQC/fastqc',
        raw_dir = IN_PATH + '/FastQC/raw/{sample}',
        clean_dir = IN_PATH + '/FastQC/clean/{sample}',
    run:
        shell('{params.FASTQC} --threads {threads} --extract -f fastq {input.raw_R1} {input.raw_R2} -o {params.raw_dir} > {log.raw} 2>&1 ')
        shell('{params.FASTQC} --threads {threads} --extract -f fastq {input.clean_R1} {input.clean_R2} -o {params.clean_dir} > {log.clean} 2>&1 ')



rule trim_QC_stats:
    input:
        raw_R1 = IN_PATH + '/raw/{sample}.R1.fastq.gz',
        raw_R2 = IN_PATH + '/raw/{sample}.R2.fastq.gz',
        clean_R1 = IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz',
        clean_R2 = IN_PATH + '/clean/{sample}.clean.paired.R2.fq.gz',
        clean = rules.fastqc.output.clean_R1,
    output:
        stat = IN_PATH + '/FastQC/{sample}_trim_stats.xls',
    log:
        IN_PATH + "/log/{sample}_trim_qc_stats.log",
    params:
        CleanDataStats = SRC_DIR + '/CleanDataStats.py',
        indir = IN_PATH + '/FastQC/clean',
    run:
        shell('python {params.CleanDataStats} --raw_R1 {input.raw_R1} --raw_R2 {input.raw_R2}  --clean_R1 {input.clean_R1}  --clean_R2  {input.clean_R2} --fastqdir {params.indir} --sample {wildcards.sample} --output {output.stat} >{log} 2>&1 ')



rule MergeQC:
    input:
        stat = expand(rules.trim_QC_stats.output.stat, sample=SAMPLES),
    output:
        stat = IN_PATH + '/FastQC/Samples_trim_stats.xls',
    run:
        INPUTS = input.stat
        print(INPUTS)
        for s in range(len(INPUTS)):
            if s == 0:
                cmd = "sed -n '1,2p' %s > %s " % (INPUTS[s], output.stat)
                os.system(cmd)
            else:
                cmd = "sed -n '2p' %s >> %s " % (INPUTS[s], output.stat)
                os.system(cmd)

