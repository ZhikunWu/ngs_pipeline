#!/usr/bin/env python
from __future__ import division
import argparse
import sys
import os
import collections

#usage: python ../SamtoolsBamStats.py --input F2_realigned_bam_stat.xls --out Samples_bam_stats.xls

#Author: Zhikun Wu
#Date: 2018.02.06
#Function: Get the statisctcs for bam file based on the result of samtools flagstat.

def bam_stats(Total, stat_file):
    in_h = open(stat_file, 'r')
    for line in in_h:
        lines = line.strip().split()
        # if 'in total' in line:
        #   Total = int(lines[0])
        if 'mapped (' in line:
            Mapped = int(lines[0])
        elif 'duplicates' in line:
            Duplicate = int(lines[0])
        # elif 'secondary' in line:
        #   Mutiple = int(lines[0])
        elif 'properly paired' in line:
            Properly = int(lines[0])
        # try:
    Mapped_ratio = Mapped / Total * 100
    Duplicate_ratio = Duplicate / Total * 100
        # Mutiple_ratio = Mutiple / Total * 100
    Properly_ratio = Properly /Total * 100
    Unmapped = Total - Mapped
    Unmapped_ratio = 100 - Mapped_ratio
    return Duplicate, Duplicate_ratio, Unmapped, Unmapped_ratio, Mapped, Mapped_ratio, Properly, Properly_ratio
        # except NameError:
        #   print('Please check whether have the mapped reads.')


def mutiple_bam_file(bam_files, clean_file, out_file):
    ### parse clean summary result
    SampleClean = collections.defaultdict()
    clean_h = open(clean_file, "r")
    headers = clean_h.readline().strip().split("\t")
    clean_index = headers.index("clean_read")
    for line in clean_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        reads = int(lines[clean_index])
        SampleClean[sample] = reads
    clean_h.close()
    ### stats of mapping reads
    out_h = open(out_file, 'w')
    out_h.write('Sample\tTotal_reads\tDuplicates(%)\tUnmapped(%)\tMapped(%)\tProperly_pair(%)\n')
    files = bam_files.split(',')
    for f in files:
        f = f.strip()
        f_base = os.path.basename(f)
        if "_realigned.bam" in f_base:
            f_base = f_base.strip("_realigned.bam")
        else:
            f_base = f_base.split("_")[0] # "_".join(f_base.split("_")[:-1])
        try:
            Total = SampleClean[f_base]
        except KeyError:
            print("Please check whether sample base %s is in the bam file %s." % (f_base, bam_files))

        Duplicate, Duplicate_ratio, Unmapped, Unmapped_ratio, Mapped, Mapped_ratio, Properly, Properly_ratio = bam_stats(Total, f)
        Duplicate_ratio = '%.3f' % Duplicate_ratio
        Unmapped_ratio = '%.3f' % Unmapped_ratio
        Mapped_ratio = '%.3f' % Mapped_ratio
        Properly_ratio = '%.3f' % Properly_ratio
        Sample = os.path.basename(f).strip('_bam_stats.xls')
        out_h.write('%s\t%d\t%d(%s)\t%d(%s)\t%d(%s)\t%d(%s)\n' % (Sample, Total, Duplicate, Duplicate_ratio, Unmapped, Unmapped_ratio, Mapped, Mapped_ratio, Properly, Properly_ratio))
    out_h.close()

def main(argv):
    parser = argparse.ArgumentParser(description='Get the statisctcs for bam file based on the result of samtools flagstat.')
    parser.add_argument('-i', '--input', help='The input files which are separated with ",".')
    parser.add_argument("-c", "--cleanSummary", help="The clean summary file.")
    parser.add_argument('-o', '--out', help='The output file.')
    args = parser.parse_args(argv)
    mutiple_bam_file(args.input, args.cleanSummary, args.out)

if __name__ == '__main__':
    main(sys.argv[1:])

