#!/usr/bin/env python 

import pandas as pd

# turn counts table into dataframe
# counts table contains amplicon name and sample name
countstable = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_counts.tsv" # open the counts table as a dataframe
countsdf = pd.read_csv(countstable, sep='\t')

# extract the Pedinophyceae amplicons from this dataframe
PedinoFile = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.uniq.txt" # open the list of Pedino amplicons
with open(PedinoFile, 'r') as file:
    PedinoList = file.read().splitlines()

PedinoDF = countsdf[countsdf['amplicon'].isin(PedinoList)] # get the Pedino subset of the counts dataframe

output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.counts.tsv"

PedinoDF.to_csv(output_file_path, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

# extract information about samples from eukbank_18S_V4_samples.tsv using the sample names from extracted from the counts table
sample_names = PedinoDF['sample'].tolist() # put all the sample names to a list

samplestable = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_samples.tsv" # open the sample table as a dataframe
samplesdf = pd.read_csv(samplestable, sep='\t')

PedinoSampleDF = samplesdf[samplesdf['sample'].isin(sample_names)] # get the Pedino subset of the sample dataframe




