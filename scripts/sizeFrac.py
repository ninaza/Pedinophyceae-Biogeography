#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:16:54 2024

@author: ninpo556
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt

def main(samples_file, output_file):
    # open the group specific sample table
    samples_table = samples_file # open the sample table as a dataframe
    samples_df = pd.read_csv(samples_table, sep='\t')
    
    size_fraction_df = samples_df[["sample", "size_fraction_lower_threshold", "size_fraction_upper_threshold"]].copy()
    size_fraction_df = size_fraction_df.set_index('sample').rename(columns={
        'size_fraction_lower_threshold': 'size_fraction_lower_threshold [µm]',
        'size_fraction_upper_threshold': 'size_fraction_upper_threshold [µm]'
    })

    size_bins = [-float('inf'), 1.9, 20, 200, float('inf')]
    size_labels = ['pico', 'nano', 'micro', 'meso']
    size_fraction_df['size_fraction_lower'] = pd.cut(size_fraction_df['size_fraction_lower_threshold [µm]'], bins=size_bins, labels=size_labels)
    size_fraction_df['size_fraction_upper'] = pd.cut(size_fraction_df['size_fraction_upper_threshold [µm]'], bins=size_bins, labels=size_labels)

    output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank_sizefractions.tsv"
    size_fraction_df.to_csv(output_file_path, sep='\t')

    # count of samples found in each category
    pico_count = size_fraction_df['size_fraction_lower'].value_counts().get('pico', 0)
    proportion_pico_water = pico_count / len(size_fraction_df) * 100
    print(proportion_pico_water)

    nano_count = size_fraction_df['size_fraction_lower'].value_counts().get('nano', 0)
    proportion_nano_water = nano_count / len(size_fraction_df) * 100
    print(proportion_nano_water)

    micro_count = size_fraction_df['size_fraction_lower'].value_counts().get('micro', 0)
    proportion_micro_water = micro_count / len(size_fraction_df) * 100
    print(proportion_micro_water)

    meso_count = size_fraction_df['size_fraction_lower'].value_counts().get('meso', 0)
    proportion_meso_water = meso_count / len(size_fraction_df) * 100
    print(proportion_meso_water)

    # barplot
    proportions = [proportion_pico_water, proportion_nano_water, proportion_micro_water, proportion_meso_water]
    categories = ['pico', 'nano', 'micro', 'meso']

    # Plot the bar plot
    colors = ['mediumaquamarine', 'cornflowerblue', 'orchid', 'orange']
    plt.figure(figsize=(12, 8), dpi=300)
    bar_width = 0.4
    plt.grid(color='white',linewidth=0.3, zorder=0)

    # Plot the bars
    plt.bar(categories, proportions, color=colors, width=bar_width, zorder=3)

    # Set background color to light grey
    plt.gca().set_facecolor('gainsboro')

    # Add labels and title
    plt.xlabel('SaSize Fraction')
    ## maybe possible to make automatic fill of target name?
    plt.ylabel('(%) Samples with target presence')
    plt.title('Proportion of Target Group in Each Size Fraction')
    
    #save plot
    output_plot_path = str(output_file + "_sizeFractions.png")
    plt.savefig(output_plot_path, bbox_inches='tight')
    
    print("Plot saved to file.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the size fraction in which your group of interest is found.')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()