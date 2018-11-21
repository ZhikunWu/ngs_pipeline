rule bowtie2:
    input:
        R1 = IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.clean.paired.R2.fq.gz',
    output:
        sam = temp(IN_PATH + '/mapping/{sample}/{sample}.sam'),
    threads:
        THREADS
    params:
        bowtie_ref = config["bowtie_ref"],
        maximum_fragment = config["maximum_fragment"],
        multimapping = 4
    log:
        IN_PATH + "/log/{sample}_bowtie2.log",
    run:
        ### -k {params.multimapping}
        shell("bowtie2  --very-sensitive -X {params.maximum_fragment} --local --mm --threads {threads} -x {params.bowtie_ref} -1 {input.R1} -2 {input.R2} > {output.sam} 2>{log}")

rule sam2bam:
    input:
        sam = IN_PATH + '/mapping/{sample}/{sample}.sam',
    output:
        bam = temp(IN_PATH + '/mapping/{sample}/{sample}.bam'),
    # params:
    #     assign_multimappers = SRC_DIR + "/assign_multimappers.py",
    run:
        ### Non-unique alignments, -q 10
        shell("samtools view -Sb -q 10 {input.sam} | samtools sort  -  > {output.bam}")


rule RemoveMitochondrial:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}.bam',
    output:
        filt_bam = temp(IN_PATH + '/mapping/{sample}/{sample}_filt_mitochrondrial.bam'),
    params:
        mitochondrial = config["mitochondrial"],
    run:
        shell("samtools index -b {input.bam}")
        shell("samtools idxstats {input.bam} | cut -f 1 | grep -v {params.mitochondrial} | xargs samtools view -b {input.bam} > {output.filt_bam}")


rule MarkDuplicates:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_filt_mitochrondrial.bam',
    output:
        bam = temp(IN_PATH + '/mapping/{sample}/{sample}_deduplicated.bam'),
        metrics = IN_PATH + '/mapping/{sample}/{sample}_dup_metrics.txt',
    params:
        PICARD = config['PICARD'],
    log:
        IN_PATH + "/log/{sample}.markDuplicate.log",
    run:
        ###  removing PCR duplicates
        shell('java -Xmx50g -jar {params.PICARD} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true >{log} 2>&1')

rule AddGroup:
    input:
        bam = rules.MarkDuplicates.output.bam,
    output:
        bam = IN_PATH + '/mapping/{sample}/{sample}_final.bam',
    params:
        PICARD = config['PICARD'],
    log:
        IN_PATH + "/log/{sample}.addGroup.log",
    run:
        shell('java -jar {params.PICARD} AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.bam} SORT_ORDER=coordinate RGID=1 RGPL=Illumina RGLB={wildcards.sample} RGPU={wildcards.sample} RGSM={wildcards.sample} >{log} 2>&1')


rule BAMIndex:
    input:
        bam = rules.AddGroup.output.bam,
    output:
        bai = IN_PATH + '/mapping/{sample}/{sample}_final.bam.bai',
    params:
        PICARD = config['PICARD'],
    log:
        IN_PATH + "/log/{sample}.bamIndex.log",
    run:
        shell('java -jar {params.PICARD} BuildBamIndex INPUT={input.bam} OUTPUT={output.bai} >{log} 2>&1')


rule BamStats:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_final.bam',
    output:
        stats = IN_PATH + '/mapping/{sample}/{sample}_bam_stats.xls',
    run:
        shell('samtools flagstat {input.bam} > {output.stats}')

rule MergeBamStats:
    input:
        stats = expand(rules.BamStats.output.stats, sample=SAMPLES),
        clean = IN_PATH + '/FastQC/Samples_trim_stats.xls',
    output:
        stats = IN_PATH + '/mapping/Samples_bam_stats.xls',
    params:
        SamtoolsBamStats = SRC_DIR + '/SamtoolsBamStats.py',
    log:
        IN_PATH + "/log/MergeBamStats.log",
    run:
        stats_files = ','.join(input.stats)
        shell('python {params.SamtoolsBamStats} --input {stats_files} --cleanSummary {input.clean} --out {output.stats} >{log} 2>&1')






