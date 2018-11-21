# rule bam2bed:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}_sorted.bam'
#     output:
#         bed = IN_PATH + '/mapping/{sample}.bed',
#     run:
#         shell("bedtools bamtobed -i {input.bam} > {output.bed}")


# ### One sample
# rule MACS:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}_final.bam',
#     output:
#         summit = IN_PATH + "/peaks/{sample}/NA_peaks.narrowPeak",
#     params:
#         outdir = IN_PATH + "/peaks/{sample}",
#         qvalue = config["qvalue"],
#         shift = config["shift"],
#         extsize = config["extsize"],
#         gsize = config["gsize"],
#     log:
#         IN_PATH + "/log/MACS_{sample}.log"
#     run:
#         ### macs2 just can install in environment python 2.7
#         shell("source activate py27 && macs2 callpeak --verbose 3 --treatment {input.bam} -g {params.gsize} -B -q {params.qvalue}  --nomodel --shift {params.shift} --extsize {params.extsize} --nolambda --keep-dup all -f BAM --outdir {params.outdir} --call-summits >{log} 2>&1")


def get_bam_file(file_list, base1, base2):
    bases = [os.path.basename(f).rstrip("_final.bam") for f in file_list]
    file1 = None
    file2 = None
    for i, j in zip(bases, file_list):
        if base1 == i:
            file1 = j
        if base2 == i:
            file2 = j
    return file1, file2


rule MACSCompare:
    input:
        bam = expand(IN_PATH + '/mapping/{sample}/{sample}_final.bam', sample=SAMPLES),
    output:
        summit = IN_PATH + "/peaks/{compare}/NA_peaks.narrowPeak",
    params:
        outdir = IN_PATH + "/peaks/{compare}",
        qvalue = config["qvalue"],
        shift = config["shift"],
        extsize = config["extsize"],
        gsize = config["gsize"],
    log:
        IN_PATH + "/log/MACS_{compare}.log"
    run:
        ### macs2 just can install in environment python 2.7
        BAMS = input.bam
        compare = wildcards.compare
        sample1 = compare.split("___")[0]
        sample2 = compare.split("___")[1]
        file1, file2 = get_bam_file(BAMS, sample1, sample2)
        shell("source activate py27 && macs2 callpeak --verbose 3 --treatment {file1} --control {file2} -g {params.gsize} -B -q {params.qvalue}  --nomodel --shift {params.shift} --extsize {params.extsize} --nolambda --keep-dup all -f BAM --outdir {params.outdir} --call-summits >{log} 2>&1")
