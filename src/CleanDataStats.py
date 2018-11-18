#!/usr/bin/env python
from __future__ import division
import gzip
import itertools
import collections
import argparse
import os

#usage: python ../CleanDataStats.py --raw_R1 raw/SRR6414595.R1.fq.gz --raw_R2 raw/SRR6414595.R2.fq.gz  --clean_R1 clean/SRR6414595.clean.paired.R1.fq.gz  --clean_R2  clean/SRR6414595.clean.paired.R2.fq.gz --fastqdir QC/clean --sample SRR6414595 --output SRR6414595._clean_stats.xls

'''
:Author: Zhikun Wu
:Email: 598466208@qq.com
:Date: 2018.03.02
:Function: Statistics for the read number and base for raw and clean data, and output the result of fastqc for clean data
'''

def fastq_read_base(fastq_file):
    """
    Get the read number and base number for the fastq file.
    """
    BaseNumber = 0
    if fastq_file.strip().endswith('gz'):
        in_h = gzip.open(fastq_file, 'r')
    else:
        in_h = open(fastq_file, 'r')
    i = 0
    lines = in_h.readlines()
    seqNumber = len(lines)
    ReadNumber = int(seqNumber/4)
    for line in itertools.islice(lines, 0, None, 4):
        seq = lines[i+1].strip()
        seqLen = len(seq)
        BaseNumber += seqLen
        i += 4
    in_h.close()
    return ReadNumber, BaseNumber

def each_trim_stat(R1_fq, R2_fq, R1_clean_paired_fq, R2_clean_paired_fq):
    '''Statistics for the read and base number for before and after trimming'''
    RawReadR1, RawBaseR1 = fastq_read_base(R1_fq)
    RawReadR2, RawBaseR2 = fastq_read_base(R2_fq)
    CleanReadR1, CleanBaseR1 = fastq_read_base(R1_clean_paired_fq)
    CleanReadR2, CleanBaseR2 = fastq_read_base(R2_clean_paired_fq)
    RawRead = RawReadR1 + RawReadR2
    CleanRead = CleanReadR1 + CleanReadR2
    try:
        Clean2RawReadRatio = CleanRead / RawRead * 100
        Clean2RawReadRatio = '%.2f' % Clean2RawReadRatio
    except KeyError:
        print('Please check whether the number of RawRead is zero.')
        sys.exit(1)
    RawBase = RawBaseR1 + RawBaseR2
    CleanBase = CleanBaseR1 + CleanBaseR2
    try:
        Clean2RawBaseRatio = CleanBase / RawBase * 100
        Clean2RawBaseRatio = '%.2f' % Clean2RawBaseRatio
    except KeyError:
        print('Please check whether the number of RawBase is zero.')        
        sys.exit(1)
    StatList = [str(RawRead), str(CleanRead), Clean2RawReadRatio, str(RawBase), str(CleanBase), Clean2RawBaseRatio]
    return StatList


def threshold_quality(q_hash,thr):
    '''Calculate the percentage of quality above the threshold for the reads.
    Usually thr1 and thr2 are 20 and 30, the headis Q20(%) and Q30(%).
    '''
    sum = 0
    above_thr = 0
    for k in q_hash.keys():
        sum += float(q_hash[k])
        if float(k) >= float(thr):
            above_thr += float(q_hash[k])
    ratio = 0
    try:
        ratio =  above_thr / sum  * 100
    except:
        print("Please Check the Per sequence quality scores result fastqc_data.txt")
    return round(ratio,1)
        
def fastqc_summary(qc_file,thr1,thr2):
    '''Summary the quality control result after performing trimmomatic and fastqc for each fastq(single) file.
        Usually thr1 and thr2 are 20 and 30. 
        used function: threshold_quality
    '''
    group_result = collections.Counter()
    quality_hash = collections.Counter()
    qc_in = open(qc_file,'r')
    for line in qc_in:
        if line.strip().startswith("%GC"):
            group_result['%GC'] = round(float(line.strip().split("\t")[1]),1)
           #Total Deduplicated Percentage^I30.876562981005602$
        elif line.strip().startswith("#Total"):
            group_result['Deduplicated_Percentage'] = round(float(line.strip().split("\t")[1]),1)
        elif line.strip().startswith(">>Per sequence quality"):
            while True:
                line = qc_in.__next__()
                if line.strip().startswith("#"):
                    continue
                elif line.strip().startswith(">>"):
                    break
                else:
                    qua_count = line.strip().split("\t")
                    quality_hash[qua_count[0]] = qua_count[1]
            threshold_quality_1 = threshold_quality(quality_hash,thr1)
            threshold_quality_2 = threshold_quality(quality_hash,thr2)
            group_result["Q"+str(thr1)+"(%)"] = threshold_quality_1
            group_result["Q"+str(thr2)+"(%)"] = threshold_quality_2
    return group_result


def two_hash_ave(h1,h2):
    '''Get the average value for two hashed those have the same keys'''
    new_h = collections.Counter()
    for i in h1.keys():
        try:
            ave = (float(h1[i]) + float(h2[i])) /2
            new_h[i] = round(ave, 1)
        except:
            print("Please check whether get the QC result from the file fastqc_data.txt!")
    return new_h
        
def fastqc_stat(sample_dir,sample_base, thr1=20, thr2=30):
    '''Running FastQC for fastq file, and return the hash contain the stats for qc result.
    '''
    sample_R1 = os.path.join(sample_dir, sample_base, sample_base + '.clean.paired.R1_fastqc/fastqc_data.txt')
    sample_R2 = os.path.join(sample_dir, sample_base, sample_base + '.clean.paired.R2_fastqc/fastqc_data.txt')
    if not (os.path.exists(sample_R1) and os.path.exists(sample_R2)):
        print("To make sure the fastqc out files %s or %s exist!" % (sample_R1,sample_R2))
    stat_R1 = fastqc_summary(sample_R1,thr1,thr2)
    stat_R2 = fastqc_summary(sample_R2,thr1,thr2)
    stat_R1_R2 = two_hash_ave(stat_R1,stat_R2)
    valueList = []
    valueList.append(str(stat_R1_R2['Q20(%)']))
    valueList.append(str(stat_R1_R2['Q30(%)']))
    valueList.append(str(stat_R1_R2['Deduplicated_Percentage']))
    valueList.append(str(stat_R1_R2['%GC']))
    return valueList

def merge_read_stat_fastqc(R1_fq, R2_fq, R1_clean_paired_fq, R2_clean_paired_fq, sample_dir,sample_base, thr1,thr2, out_file):
    StatList = each_trim_stat(R1_fq, R2_fq, R1_clean_paired_fq, R2_clean_paired_fq)
    valueList = fastqc_stat(sample_dir,sample_base, thr1,thr2)
    out_h = open(out_file, 'w')
    out_h.write("Sample\tQ20(%)\tQ30(%)\tDeduplicated_Percentage\tGC(%)\traw_read\tclean_read\tclean2raw_read_ratio(%)\traw_base\tclean_base\tclean2raw_base_ratio(%)\n")
    out_h.write('%s\t%s\t%s\n' % (sample_base, '\t'.join(valueList), '\t'.join(StatList)))
    out_h.close()


def main():
    formatter_class = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="summarize the results of fastq files", formatter_class=formatter_class)
    parser.add_argument('-d','--fastqdir',help='The dir for the clean fastqc out directory.')
    # parser.add_argument('-i','--input',help='The input file which is string seperated with comma for many fastqc_data.txt ')
    parser.add_argument('-o','--output',help='The statistics of fastqc results.')
    parser.add_argument('-f','--raw_R1',help='The list of read and base statistics files which were seperated with ')
    parser.add_argument('-r','--raw_R2',help='The list of read and base statistics files which were seperated with ')
    parser.add_argument('-c','--clean_R1',help='The list of read and base statistics files which were seperated with ')
    parser.add_argument('-e','--clean_R2',help='The list of read and base statistics files which were seperated with ')
    parser.add_argument('-s', '--sample', help='The sample base name.')
    args = parser.parse_args()
    merge_read_stat_fastqc(args.raw_R1, args.raw_R2, args.clean_R1, args.clean_R2, args.fastqdir, args.sample, 20, 30, args.output)

if __name__ == '__main__':
    main()
