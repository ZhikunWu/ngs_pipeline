#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import sys
import yaml


IN_PATH = config['project_path']
THREADS = config["THREADS"]
SAMPLES = config['SAMPLES']
BIN_DIR = config["BIN_DIR"]
PIP_DIR = config["PIP_DIR"]
RULE_DIR = PIP_DIR + "/../rule"
SRC_DIR = PIP_DIR + "/../src"

include: RULE_DIR + "/QualityControl.rule.py"
include: RULE_DIR + "/ReadsMapping.rule.py"


rule all:
    input:
        ################### Quality Control #################################
        expand(IN_PATH + '/clean/{sample}.clean.paired.R1.fq.gz', sample=SAMPLES),
        expand(IN_PATH + '/FastQC/raw/{sample}/{sample}.R1_fastqc.zip', sample=SAMPLES),
        expand(IN_PATH + '/FastQC/{sample}_trim_stats.xls', sample=SAMPLES),
        IN_PATH + '/FastQC/Samples_trim_stats.xls',
        ########################################################################
        ################################# mapping ##############################
        expand(IN_PATH + '/mapping/{sample}.sam', sample=SAMPLES),
        expand(IN_PATH + '/mapping/{sample}.sorted.bam', sample=SAMPLES),
        expand(IN_PATH + '/mapping/{sample}.bed', sample=SAMPLES),
