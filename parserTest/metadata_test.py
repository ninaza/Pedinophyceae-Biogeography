#!/usr/bin/env python3

import pandas as pd
import argparse

## usage: python metadata.py -g group_file.txt -c counts_file.tsv -s samples_file.csv -o output_file ##

def main(group_file, counts_file, samples_file, output_file):
    # turn counts table into dataframe
    # counts table contains amplicon name and sample name
    counts_table = counts_file # open the counts table as a dataframe
    counts_df = pd.read_csv(counts_table, sep='\t')

    # extract the Pedinophyceae amplicons from this dataframe
    group_file = group_file # open the list of Pedino amplicons
    with open(group_file, 'r') as file:
        group_list = file.read().splitlines()

    group_df = counts_df[counts_df['amplicon'].isin(group_list)] # get the Pedino subset of the counts dataframe

    output_file_path = str(output_file + "_counts.tsv")

    group_df.to_csv(output_file_path, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

    # extract information about samples from eukbank_18S_V4_samples.tsv using the sample names from extracted from the counts table
    sample_names = group_df['sample'].tolist() # put all the sample names to a list

    samples_table = samples_file # open the sample table as a dataframe
    samples_df = pd.read_csv(samples_table, sep='\t')

    group_sample_df = samples_df[samples_df['sample'].isin(sample_names)] # get the Pedino subset of the sample dataframe

    output_file_path2 = str(output_file + "_samples.tsv")

    group_sample_df.to_csv(output_file_path2, sep='\t', index=False) # save the Pedino counts dataframe as tsv file

    print("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subset the count and sample table for your desired group.')
    parser.add_argument('-g', '--group', type=str, help='Path to group file')
    parser.add_argument('-c', '--counts', type=str, help='Path to counts file')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.group, args.counts, args.samples, args.output)