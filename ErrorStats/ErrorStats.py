#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:18:47 2019

@author: miaaltieri

ran with both:
    deepvariant_hg002_grch37.vcf.gz
    friday_hg002_grch37.vcf.gz 
"""
from pysam import VariantFile
import seaborn as sns
import pandas as pd
import sys

# Step 0: Check inputs
if len(sys.argv) != 2:
    print("EROR: missing the required vcf file")
    print("      file should be: of the form .vcf.gz")
    sys.exit()
    
filename = sys.argv[1]
# End Step 0

# Step 1: Split the VCF into two files 
vcf_in = VariantFile(filename)  # auto-detect input format
vcf_no_calls = VariantFile("NoCalls.vcf", 'w', header=vcf_in.header)
vcf_genotype = VariantFile("Genotype.vcf", 'w', header=vcf_in.header)


for rec in vcf_in.fetch():
    if rec.samples['TRUTH']['BD'] == 'FN':
        if rec.samples['QUERY']['BLT'] == 'nocall':
            vcf_no_calls.write(rec)
        else:
            vcf_genotype.write(rec)
            
# End Step 1


# Step 2: Iterate through nocalls and genotype and create dictionaries to
# count which chromosomes have which errors
vcf_no_calls = VariantFile("NoCalls.vcf")
vcf_genotype = VariantFile("Genotype.vcf")

no_calls_errors = {}
genotype_errors = {}

for rec in vcf_no_calls.fetch():
    chrom = int(rec.chrom)
    if chrom not in no_calls_errors:
        no_calls_errors[chrom] = 0
    no_calls_errors[chrom] = no_calls_errors[chrom] + 1
    
for rec in vcf_genotype.fetch():
    chrom = int(rec.chrom)
    if chrom not in genotype_errors:
       genotype_errors[chrom] = 0
    genotype_errors[chrom] =  genotype_errors[chrom] + 1
    
no_call_total = 0
for key in no_calls_errors:
    no_call_total += no_calls_errors[key]
    
genotype_total = 0
for key in genotype_errors:
    genotype_total += genotype_errors[key]
    
# End Step 2

# Step 3: create seaborn graphs based on these results
# VERSION A: This gives you two different graphs for each type
data_as_dict = {}
chrom_vals = []
error_vals = []


for i in range (1,23):
    chrom_vals.append(i)
    error_vals.append(no_calls_errors[i])

data_as_dict['error count']=error_vals
data_as_dict['chromosome']=chrom_vals
    
df = pd.DataFrame(data_as_dict)

sns.set(font_scale=1.5)
g = sns.catplot(x="chromosome", y="error count", kind="bar", data=df, ci=.95);
g.fig.set_size_inches(14,8)
g.set(ylim=(0,600))
axes = g.axes.flatten()
axes[0].set_title('Error Rates for No Calls')


# same thing but for the genotype errors
data_as_dict = {}
chrom_vals = []
error_vals = []


for i in range (1,23):
    chrom_vals.append(i)
    error_vals.append(genotype_errors[i])

data_as_dict['error count']=error_vals
data_as_dict['chromosome']=chrom_vals
    
df = pd.DataFrame(data_as_dict)

sns.set(font_scale=1.5)
g = sns.catplot(x="chromosome", y="error count", kind="bar", data=df, ci=.95);
g.fig.set_size_inches(14,8)
g.set(ylim=(0,600))
axes = g.axes.flatten()
axes[0].set_title('Error Rates for Genotypes')

# VERSION B: This gives you one graph and compares the two together
data_as_dict = {}
chrom_vals = []
error_vals = []
error_type = []


for i in range (1,23):
    chrom_vals.append(i)
    chrom_vals.append(i)
    error_vals.append(no_calls_errors[i])
    error_vals.append(genotype_errors[i])
    error_type.append('no calls')
    error_type.append('genotype')
    
data_as_dict['chromosome']=chrom_vals
data_as_dict['error count']=error_vals
data_as_dict['error type']=error_type

    
df = pd.DataFrame(data_as_dict)

sns.set(font_scale=1.5)
g = sns.catplot(x="chromosome", y="error count", hue='error type', kind="bar", data=df, ci=.95);
g.fig.set_size_inches(20,8)
g.set(ylim=(0,600))
axes = g.axes.flatten()
axes[0].set_title('Comparing No Call with Genotype Error Rates')

# VERSION C: This gives you the same thing as B, but instead it is as a percentage
# of the errors made

# VERSION B: This gives you one graph and compares the two together
data_as_dict = {}
chrom_vals = []
error_vals = []
error_type = []


for i in range (1,23):
    chrom_vals.append(i)
    chrom_vals.append(i)
    error_vals.append(no_calls_errors[i]/no_call_total*100)
    error_vals.append(genotype_errors[i]/genotype_total*100)
    error_type.append('no calls')
    error_type.append('genotype')
    
data_as_dict['chromosome']=chrom_vals
data_as_dict['error percent']=error_vals
data_as_dict['error type']=error_type

    
df = pd.DataFrame(data_as_dict)

sns.set(font_scale=1.5)
g = sns.catplot(x="chromosome", y="error percent", hue='error type', kind="bar", data=df, ci=.95);
g.fig.set_size_inches(20,8)
g.set(ylim=(0,15))
axes = g.axes.flatten()
axes[0].set_title('Comparing No Call with Genotype Error Rates by Percent of Total Errors')
# End Step 3
