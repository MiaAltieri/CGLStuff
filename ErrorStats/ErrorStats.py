#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:18:47 2019

@author: miaaltieri
"""
from pysam import VariantFile

# Step 1: Split the VCF into two files 

vcf_in = VariantFile("friday_hg002_grch37.vcf.gz")  # auto-detect input format
vcf_no_calls = VariantFile("NoCalls.vcf", 'w', header=vcf_in.header)
vcf_genotype = VariantFile("Genotype.vcf", 'w', header=vcf_in.header)

i=0
for rec in vcf_in.fetch():
    if i == 0:
        continue
    if i == 2:
        break
    vcf_no_calls.write(rec)
    vcf_genotype.write(rec)
    
    i+=1