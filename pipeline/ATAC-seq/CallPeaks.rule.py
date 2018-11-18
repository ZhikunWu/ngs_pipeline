# rule bam2bed:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}_sorted.bam'
#     output:
#         bed = IN_PATH + '/mapping/{sample}.bed',
#     run:
#         shell("bedtools bamtobed -i {input.bam} > {output.bed}")

rule MACS:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_final.bam',
    output:
        summit = IN_PATH + "/peaks/{sample}/NA_summits.bed",
    params:
        outdir = IN_PATH + "/peaks/{sample}",
        qvalue = config["qvalue"],
        shift = config["shift"],
        extsize = config["extsize"],
        gsize = config["gsize"],
    log:
        IN_PATH + "/log/MACS_{sample}.log"
    run:
        ### macs2 just can install in environment python 2.7
        shell("source activate py27 && macs2 callpeak --verbose 3 --treatment {input.bam} -g {params.gsize} -B -q {params.qvalue}  --nomodel --shift {params.shift} --extsize {params.extsize} --nolambda --keep-dup all -f BAM --outdir {params.outdir} --call-summits >{log} 2>&1")
