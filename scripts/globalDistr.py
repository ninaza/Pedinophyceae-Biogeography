#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:22:14 2024

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
    
    lon_lat_all_df = samples_df[["sample", "latitude", "longitude"]].copy()

    lon_lat_all_df['state'] = lon_lat_all_df['sample'].isin(sample_names).map({True: 'yes', False: 'no'})

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

    # Plot map on axis
    countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    countries.plot(color="lightgrey", ax=ax)

    # Plot points with state "no" as crosses
    no_state_df = lon_lat_all_df[lon_lat_all_df['state'] == 'no']
    ax.scatter(no_state_df['longitude'], no_state_df['latitude'], marker='x', linewidth=0.4, s=10, color='dimgrey', label='No')

    # Plot points with state "yes" as dots with size corresponding to percentage
    yes_state_df = lon_lat_all_df[lon_lat_all_df['state'] == 'yes']
    ax.scatter(yes_state_df['longitude'], yes_state_df['latitude'], marker='o', color='cornflowerblue', label='Yes', alpha=0.5)

    # Set labels and title
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Add legend
    custom_handles = [ax.scatter([], [], marker='o', color='cornflowerblue', s=30, label='Present'),
                      ax.scatter([], [], marker='x', color='dimgrey', s=30, linewidth=0.4, label='Absent')]

    ax.legend(handles=custom_handles, loc='lower left', bbox_to_anchor=(1, 0.6))

    #save plot
    output_plot_path = str(output_file + "_globalDistribution.png")
    plt.savefig(output_plot_path, bbox_inches='tight')
    
    print("Plot saved to file.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the global distribution of your group of interest.')
    parser.add_argument('-c', '--counts', type=str, help='Path to counts file')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.counts, args.samples, args.output)