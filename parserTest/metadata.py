#!/usr/bin/env python 

import pandas as pd

# turn counts table into dataframe
# counts table contains amplicon name and sample name
counts_table = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_counts.tsv" # open the counts table as a dataframe
counts_df = pd.read_csv(counts_table, sep='\t')

# extract the Pedinophyceae amplicons from this dataframe
group_file = "Pedino-EukBank.uniq.txt" # open the list of Pedino amplicons
with open(group_file, 'r') as file:
    group_list = file.read().splitlines()

group_df = counts_df[counts_df['amplicon'].isin(group_list)] # get the Pedino subset of the counts dataframe

output_file1 = "Pedino_counts.tsv"

group_df.to_csv(output_file1, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

# extract information about samples from eukbank_18S_V4_samples.tsv using the sample names from extracted from the counts table
sample_names = group_df['sample'].tolist() # put all the sample names to a list

samples_table = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_samples.tsv" # open the sample table as a dataframe
samples_df = pd.read_csv(samples_table, sep='\t')

group_sample_df = samples_df[samples_df['sample'].isin(sample_names)] # get the Pedino subset of the sample dataframe

output_file2 = "Pedino_samples.tsv"
group_sample_df.to_csv(output_file2, sep='\t', index=False) # save the Pedino counts dataframe as tsv file
