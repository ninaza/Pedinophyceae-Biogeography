#!/usr/bin/env python 

import pandas as pd

# turn counts table into dataframe
# counts table contains amplicon name and sample name
counts_table = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_counts.tsv" # open the counts table as a dataframe
counts_df = pd.read_csv(counts_table, sep='\t')

# extract the Pedinophyceae amplicons from this dataframe
group_file = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.uniq.txt" # open the list of Pedino amplicons
with open(group_file, 'r') as file:
    group_list = file.read().splitlines()

group_df = counts_df[counts_df['amplicon'].isin(group_list)] # get the Pedino subset of the counts dataframe

output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.counts.tsv"

group_df.to_csv(output_file_path, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

# extract information about samples from eukbank_18S_V4_samples.tsv using the sample names from extracted from the counts table
sample_names = group_df['sample'].tolist() # put all the sample names to a list

samples_table = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_samples.tsv" # open the sample table as a dataframe
samples_df = pd.read_csv(samples_table, sep='\t')

group_sample_df = samples_df[samples_df['sample'].isin(sample_names)] # get the Pedino subset of the sample dataframe

output_file_path2 = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.tsv"

group_sample_df.to_csv(output_file_path2, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

#if __name__ == '__main__':
#    parser = ArgumentParser(prog='EukBank Metadata')
#    parser.add_argument('-c', dest='counts_table', help='counts_table', nargs=1)
#    parser.add_argument('-g', dest='group_file', help='group_file', nargs=1)
#    parser.add_argument('-s', dest='samples_table', help='samples_table', nargs=1)