#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:18:47 2019

@author: miaaltieri
"""
from pysam import VariantFile

# Step 1: Split the VCF into two files 

vcf_in = VariantFile("pfda_HG001_GRCh37_02152019_182915.vcf.gz")  # auto-detect input format
vcf_no_calls = VariantFile("NoCalls.vcf", 'w', header=vcf_in.header)
vcf_genotype = VariantFile("NoCalls.vcf", 'w', header=vcf_in.header)

for rec in vcf_in.fetch('chr1', 100000, 200000):
    vcf_no_calls.write(rec)
    vcf_genotype.write(rec)