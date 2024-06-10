#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 11:45:47 2024

@author: ninpo556
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import geopandas as gpd

def main(counts_file, samples_file, output_file):
    # open group specific counts table
    counts_table = counts_file # open the counts table as a dataframe
    counts_df = pd.read_csv(counts_table, sep='\t')
    
    # extract all sample names
    sample_names = counts_df['sample'].tolist() # put all the sample names to a list
    
    # open the general EukBank sample table
    samples_table = samples_file # open the sample table as a dataframe
    samples_df = pd.read_csv(samples_table, sep='\t')
    
    # subset sample_df to sample, latitude and longitude comlumns
    lon_lat_all_df = samples_df[["sample", "latitude", "longitude"]].copy()
    # add new column with presence and absence information
    lon_lat_all_df['state'] = lon_lat_all_df['sample'].isin(sample_names).map({True: 'yes', False: 'no'})
    
    ## calculation of relative abundance of Pedinophyceae endosymbionts per sample
    # Group by 'Sample' and calculate the mean number of reads per sample
    average_reads_per_sample = counts_df.groupby('sample')['nreads'].mean()
    average_reads_per_sample = average_reads_per_sample.rename('nreads-Pedino')

    # take sample and nreads columns from sample_df, set sample name as index and rename nreads to nreads total
    total_sample_abundance_df = samples_df[["sample", "nreads"]].copy()
    total_sample_abundance_df = total_sample_abundance_df.set_index('sample')
    total_sample_abundance_df = total_sample_abundance_df.rename(columns={'nreads': 'nreads-total'})

    # merge df with mean abundance of group per sample and df with total number of reads in sample and calulcate ratio and percentage
    merged_df = total_sample_abundance_df.merge(average_reads_per_sample, left_on='sample', right_index=True)
    merged_df['relative_abundance'] = merged_df['nreads-Pedino'] / merged_df['nreads-total']
    merged_df['abundance_percent'] = merged_df['relative_abundance'] * 100

    percentage_df = merged_df[["abundance_percent"]].copy()
    lon_lat_all_df = lon_lat_all_df.set_index('sample')
    lon_lat_all_percentage_df = pd.merge(lon_lat_all_df, percentage_df, on='sample', how='left')
    lon_lat_all_percentage_df['abundance_percent'].fillna(0, inplace=True)

    # Find the sample with the highest abundance
    max_abundance_sample = merged_df.loc[merged_df['relative_abundance'].idxmax()]
    print("Sample with the highest abundance:")
    print(max_abundance_sample)
    
    # Prompt the user for the max abundance value
    max_abundance = float(input("Please enter the max abundance value: "))

    # Calculate percentage values for the legend
    percentage_values = [max_abundance, max_abundance / 4, (max_abundance / 4) / 2.5, ((max_abundance / 4) / 2.5) / 4, (((max_abundance / 4) / 2.5) / 4) / 2.5]

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

    # Plot map on axis
    countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    countries.plot(color="lightgrey", ax=ax)

    # Plot points with state "no" as crosses
    no_state_df = lon_lat_all_percentage_df[lon_lat_all_percentage_df['state'] == 'no']
    ax.scatter(no_state_df['longitude'], no_state_df['latitude'], marker='x', linewidth=0.4, s=10, color='dimgrey', label='No')

    # Plot points with state "yes" as dots with size corresponding to percentage
    yes_state_df = lon_lat_all_percentage_df[lon_lat_all_percentage_df['state'] == 'yes']
    ax.scatter(yes_state_df['longitude'], yes_state_df['latitude'], s=yes_state_df['abundance_percent'] * 100, color='cornflowerblue', label='Yes', alpha=0.5)

    # Set labels and title
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Create legend scatter plots for reference sizes
    legend_handles = []
    for pct in percentage_values:
        size = pct * 100  # Scaling the size for better visualization
        handle = ax.scatter([], [], s=size, label=f'{pct:.2f}%', c='cornflowerblue', alpha=0.5)
        legend_handles.append(handle)

    # Add legend
    ax.legend(handles=legend_handles, title='Relative abundance', loc='lower left', bbox_to_anchor=(1, 0.5),)
    
    #save plot
    output_plot_path = str(output_file + "_averageAbundance.png")
    plt.savefig(output_plot_path, bbox_inches='tight')
    
    print("Plot saved to file.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the average abundance of your group of interest.')
    parser.add_argument('-c', '--counts', type=str, help='Path to counts file')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.counts, args.samples, args.output)
    