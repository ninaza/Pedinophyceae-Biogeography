#!/usr/bin/env python 

# plotting of the sampling spot

import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd

sample_position_df = pd.read_csv("/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.tsv", 
                              sep='\t', usecols=["latitude", "longitude"])
# initialize an axis
fig, ax = plt.subplots(figsize=(8,6))

# plot map on axis
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
countries.plot(color="lightgrey", ax=ax)

# plot points
sample_position_df.plot(x="longitude", y="latitude", kind="scatter", ax=ax)

plt.savefig('/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.png', 
            bbox_inches='tight')

# calculation of relative abundance of Pedinophyceae endosymbionts per sample

counts_table = "/Users/ninpo556/Desktop/DB/EukBank/eukbank_18S_V4_counts.tsv" # open the counts table as a dataframe
counts_df = pd.read_csv(counts_table, sep='\t')
# Group by 'Sample' and calculate the mean number of reads per sample
average_reads_per_sample = counts_df.groupby('sample')['nreads'].mean()
average_reads_per_sample = average_reads_per_sample.rename('nreads-Pedino')

# Print the resulting Series containing the average reads per sample
print("Average Reads per Sample:")
print(average_reads_per_sample)

total_sample_abundance_df = pd.read_csv("/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.tsv",
                              sep='\t', usecols=["sample", "nreads"])
total_sample_abundance_df = total_sample_abundance_df.set_index('sample')
total_sample_abundance_df = total_sample_abundance_df.rename(columns={'nreads': 'nreads-total'})

merged_df = total_sample_abundance_df.merge(average_reads_per_sample, left_on='sample', right_index=True)

merged_df['relative_abundance'] = merged_df['nreads-Pedino'] / merged_df['nreads-total']

max_division_row = merged_df.loc[merged_df['relative_abundance'].idxmax()]

output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank_nreads.tsv"

# Save the DataFrame to a TSV file
merged_df.to_csv(output_file_path, sep='\t')
