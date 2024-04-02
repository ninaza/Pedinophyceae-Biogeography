#!/usr/bin/env python 

####### DISTRIBUTION AND ABUNDANCE ############################################
# plotting of the sampling spot

import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
#import cartopy.crs as ccrs
import argparse

def main(samples_file,counts_file, output_file):
    
    #### PLOT SAMPLE LOCATIONS
    
    ## open group sample table
    samples_table = samples_file # open the sample table as a dataframe
    samples_df = pd.read_csv(samples_table, sep='\t')

    ## extract longitude and latitude data
    lon_lat_df = samples_df[["latitude", "longitude"]].copy()
    # initialize an axis
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

    # plot map on axis
    countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    countries.plot(color="lightgrey", ax=ax)

    # plot points
    lon_lat_df.plot(x="longitude", y="latitude", kind="scatter", ax=ax)
    
    output_plt1 = str(output_file + "_samples.png")
    
    plt.savefig(output_plt1, bbox_inches='tight')

    print("1) Plotting sample locations.")

    #### CALCULATE RELATIVE ABUNDANCE AND INCLUDE IN PLOT
    
    # calculation of relative abundance of Pedinophyceae endosymbionts per sample
    counts_table = counts_file # open the counts table as a dataframe
    counts_df = pd.read_csv(counts_table, sep='\t')
    # Group by 'Sample' and calculate the mean number of reads per sample
    average_reads_per_sample = counts_df.groupby('sample')['nreads'].mean()
    average_reads_per_sample = average_reads_per_sample.rename(str("nreads-" + output_file))
    
    total_sample_abundance_df = samples_df[["sample", "nreads"]].copy()
    total_sample_abundance_df = total_sample_abundance_df.set_index('sample')
    total_sample_abundance_df = total_sample_abundance_df.rename(columns={'nreads': 'nreads-total'})

    merged_df = total_sample_abundance_df.merge(average_reads_per_sample, left_on='sample', right_index=True)
    merged_df['relative_abundance'] = merged_df[str("nreads-" + output_file)] / merged_df['nreads-total']
    merged_df['abundance_percent'] = merged_df['relative_abundance'] * 100

    output_file1 = str(output_file + "_nreads.tsv")
    merged_df.to_csv(output_file1, sep='\t', index=False)
    
    print("2) Saving relative abundances to file.")
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subset the count and sample table for your desired group.')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-c', '--counts', type=str, help='Path to counts file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.samples, args.counts, args.output)
