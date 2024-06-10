#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:08:26 2024

@author: ninpo556
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import geopandas as gpd

def main(samples_file, output_file):
    # open group sample table
    sample_table = samples_file # open the counts table as a dataframe
    sample_df = pd.read_csv(sample_table, sep='\t')
    
    # assign depth lable: sunlit (0-200 m) and dark ocean (>200 m)
    depth_df = sample_df[["sample", "depth"]].copy()
    
    # get maximal depth
    max_depth = depth_df.loc[depth_df['depth'].idxmax()]
    print(max_depth)

    # Categorize the depth into two categories
    depth_df['ocean_layer'] = pd.cut(depth_df['depth'], bins=[-float('inf'), 200, float('inf')], labels=['sunlit ocean', 'dark ocean'])
    depth_df['ocean_layer'] = depth_df['ocean_layer'].cat.add_categories(['unknown'])
    depth_df['ocean_layer'].fillna('unknown', inplace=True)

    # count of samples found in sunlit ocean
    sunlit_ocean_count = depth_df['ocean_layer'].value_counts().get('sunlit ocean', 0)
    proportion_sunlit_ocean = sunlit_ocean_count / len(depth_df) * 100
    print(proportion_sunlit_ocean)

    # count of samples found in dark ocean
    dark_ocean_count = depth_df['ocean_layer'].value_counts().get('dark ocean', 0)
    proportion_dark_ocean = dark_ocean_count / len(depth_df) * 100
    print(proportion_dark_ocean)

    # plot category of each sample
    sample_loc_df = sample_df[["sample", "latitude", "longitude"]].copy()
    merged_df = pd.merge(sample_loc_df, depth_df, on="sample", how="left")

    # Initialize an axis
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

    # Plot map on axis
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world.plot(ax=ax, color='lightgrey')

    # Plot points with color based on depth
    colors = {'sunlit ocean': 'dodgerblue', 'dark ocean': 'navy', 'unknown': 'darkgrey'}  # Mapping of depth to colors
    for ocean_layer, group in merged_df.groupby('ocean_layer'):
        group.plot(kind='scatter', x='longitude', y='latitude', color=colors[ocean_layer], alpha=0.5, label=ocean_layer, ax=ax)

    # Add legend
    legend_labels = ['sunlit ocean', 'dark ocean', 'unknown']
    handles = [plt.Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=colors[ocean_layer]) for ocean_layer in legend_labels]
    ax.legend(handles, legend_labels)
    
    #save plot
    output_plot_path = str(output_file + "_depth.png")
    plt.savefig(output_plot_path, bbox_inches='tight')
    
    print("Plot saved to file.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the depth distribution of your group of interest.')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.samples, args.output)