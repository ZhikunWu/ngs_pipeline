rule bowtie2:
    input:
        R1 = IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.clean.paired.R2.fq.gz',
    output:
        sam = IN_PATH + '/mapping/{sample}.sam',
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
        shell("bowtie2  -X {params.maximum_fragment} --local --mm --threads {threads} -x {params.bowtie_ref} -1 {input.R1} -2 {input.R2} > {output.sam} 2>{log}")

rule sam2bam:
    input:
        sam = IN_PATH + '/mapping/{sample}.sam',
    output:
        # bam = IN_PATH + '/mapping/{sample}.bam',
        # filtMultip_bam = IN_PATH + '/mapping/{sample}.filtMultip.bam',
        # filt_bam = IN_PATH + '/mapping/{sample}.filt.bam',
        # fixmate = IN_PATH + '/mapping/{sample}.fixmate.tmp',
        sorted_bam = IN_PATH + '/mapping/{sample}.sorted.bam',
    params:
        assign_multimappers = SRC_DIR + "/assign_multimappers.py",
        # multimapping = 4
    run:
        shell("samtools view -Sb {input.sam} | samtools sort -  > {output.sorted_bam}")
        # shell("samtools view -F 524 -f 2 -u  {output.bam} | samtools sort -n - > {output.filtMultip_bam}")
        # shell("samtools view -h {output.filtMultip_bam} | {params.assign_multimappers} -k {params.multimapping} -paired-end | samtools view -bS - > {output.filt_bam}")
        # shell("samtools fixmate -r {output.filt_bam} {output.fixmate}")
        # shll("samtools view -F 1804 -f 2 -u {output.fixmate} | samtools sort - > {output.sorted_bam}")


rule bam2bed:
    input:
        bam = IN_PATH + '/mapping/{sample}.sorted.bam',
    output:
        bed = IN_PATH + '/mapping/{sample}.bed',
    run:
        shell("bedtools bamtobed -i {input.bam} > {output.bed}")