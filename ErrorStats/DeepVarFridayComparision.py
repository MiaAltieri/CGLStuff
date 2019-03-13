#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:03:20 2019

ran with:
    deepvariant_hg002_grch37.vcf.gz
    friday_hg002_grch37.vcf.gz 

@author: miaaltieri
"""

from pysam import VariantFile
import seaborn as sns
import pandas as pd
import sys

# Step 0: Check inputs
if len(sys.argv) != 3:
    print("EROR: missing the two required files")
    print("      first file should be: deep_variant.vcf.gz")
    print("      second file should be: friday.vcf.gz")
    sys.exit()
    
deep_variant_filename = sys.argv[1]
friday_filename = sys.argv[2]
# End Step 0

# Step 1: Split the two VCFs into three files
vcf_deep_var = VariantFile(deep_variant_filename) 
vcf_friday = VariantFile(friday_filename)
vcf_deep_var_errors = VariantFile("DeepVariantOnly.vcf", 'w', header=vcf_deep_var.header)
vcf_friday_errors = VariantFile("FridayOnly.vcf", 'w', header=vcf_friday.header)
vcf_both_errors = VariantFile("Both.vcf", 'w', header=vcf_friday.header)

deep_var = {}
friday = {}
both = {}

for rec in vcf_deep_var.fetch():
    if rec.samples['TRUTH']['BD'] == 'FN':
        # unique key is the chromosome and position within chromosome
        key = str(rec.chrom) + ' ' + str(rec.pos)
        deep_var[key] = rec

for rec in vcf_friday.fetch():
    if rec.samples['TRUTH']['BD'] == 'FN':
        # unique key is the chromosome and position within chromosome
        key = str(rec.chrom) + ' ' + str(rec.pos)
        if key in deep_var:
            del deep_var[key]
            both[key] = rec
        else: 
            friday[key] = rec


for key, val in deep_var.items():
    vcf_deep_var_errors.write(val)

for key, val in friday.items():
    vcf_friday_errors.write(val)

for key, val in both.items():
    vcf_both_errors.write(val)
    
# End Step 1

# Step 2: Iterate through nocalls and genotype and create dictionaries to
# count which chromosomes have which errors
vcf_deep_var_errors = VariantFile("DeepVariantOnly.vcf")
vcf_friday_errors = VariantFile("FridayOnly.vcf")
vcf_both_errors = VariantFile("Both.vcf")

deep_var_errors = {}
friday_errors = {}
both_errors = {}

for rec in vcf_deep_var_errors.fetch():
    chrom = int(rec.chrom)
    if chrom not in deep_var_errors:
        deep_var_errors[chrom] = 0
    deep_var_errors[chrom] = deep_var_errors[chrom] + 1
    
for rec in vcf_friday_errors.fetch():
    chrom = int(rec.chrom)
    if chrom not in friday_errors:
       friday_errors[chrom] = 0
    friday_errors[chrom] =  friday_errors[chrom] + 1
    
for rec in vcf_both_errors.fetch():
    chrom = int(rec.chrom)
    if chrom not in both_errors:
       both_errors[chrom] = 0
    both_errors[chrom] =  both_errors[chrom] + 1
# End Step 2
    
# Step 3: create seaborn graphs based on these results
# This gives you one graph and compares the three vcfs together
data_as_dict = {}
chrom_vals = []
error_vals = []
error_type = []

for i in range (1,23):
    chrom_vals.append(i)
    chrom_vals.append(i)
    chrom_vals.append(i)
    error_vals.append(deep_var_errors[i])
    error_vals.append(friday_errors[i])
    error_vals.append(both_errors[i])
    error_type.append('Deep Variant')
    error_type.append('Friday')
    error_type.append('Both')
    
data_as_dict['chromosome']=chrom_vals
data_as_dict['error count']=error_vals
data_as_dict['error type']=error_type

df = pd.DataFrame(data_as_dict)

sns.set(font_scale=1.5)
g = sns.catplot(x="chromosome", y="error count", hue='error type', kind="bar", data=df, ci=.95);
g.fig.set_size_inches(20,8)
g.set(ylim=(0,600))
axes = g.axes.flatten()
axes[0].set_title('Comparing Errors between DeepVar and Friday Models')
# End Step 3