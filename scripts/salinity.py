#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:35:51 2024

@author: ninpo556
"""

import pandas as pd
import argparse
import matplotlib.pyplot as plt

def main(samples_file, output_file):
    # open the group specific sample table
    samples_table = samples_file # open the sample table as a dataframe
    samples_df = pd.read_csv(samples_table, sep='\t')
    
    # assign salinity lable: briny (>50 ppt), saline water (30-50 ppt), brackish water (0.5-30 ppt), fresh water (0-0.5 ppt)
    salinity_df = samples_df[["sample", "salinity"]].copy()
    salinity_df.set_index('sample', inplace=True)

    # get maximal depth
    max_salinity = salinity_df.loc[salinity_df['salinity'].idxmax()]
    print(max_salinity)

    # Categorize the salinity
    bins = [-float('inf'), 0.5, 30, 50, float('inf')]
    labels = ['fresh water', 'brackish water', 'saline water', 'briny water']
    salinity_df['salinty_range'] = pd.cut(salinity_df['salinity'], bins=bins, labels=labels, right=False)
    salinity_df['salinty_range'] = salinity_df['salinty_range'].cat.add_categories(['unknown'])
    salinity_df['salinty_range'].fillna('unknown', inplace=True)

    output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank_salinity.tsv"
    salinity_df.to_csv(output_file_path, sep='\t')

    # count of samples found in each category
    briny_water_count = salinity_df['salinty_range'].value_counts().get('briny water', 0)
    proportion_briny_water = briny_water_count / len(salinity_df) * 100
    print(proportion_briny_water)

    saline_water_count = salinity_df['salinty_range'].value_counts().get('saline water', 0)
    proportion_saline_water = saline_water_count / len(salinity_df) * 100
    print(proportion_saline_water)

    brackish_water_count = salinity_df['salinty_range'].value_counts().get('brackish water', 0)
    proportion_brackish_water = brackish_water_count / len(salinity_df) * 100
    print(proportion_brackish_water)

    fresh_water_count = salinity_df['salinty_range'].value_counts().get('fresh water', 0)
    proportion_fresh_water = fresh_water_count / len(salinity_df) * 100
    print(proportion_fresh_water)

    # barplot
    proportions = [proportion_fresh_water, proportion_brackish_water, proportion_saline_water, proportion_briny_water]
    categories = ['Fresh Water', 'Brackish Water', 'Saline Water', 'Briny Water']

    # Plot the bar plot
    colors = ['lightblue', 'skyblue', 'deepskyblue', 'dodgerblue']
    plt.figure(figsize=(12, 8), dpi=300)
    bar_width = 0.4
    plt.grid(color='white',linewidth=0.3, zorder=0)

    # Plot the bars
    plt.bar(categories, proportions, color=colors, width=bar_width, zorder=3)

    # Set background color to light grey
    plt.gca().set_facecolor('gainsboro')

    # Add labels and title
    plt.xlabel('Salinity Category')
    ## maybe possible to make automatic fill of target name?
    plt.ylabel('(%) Samples with target presence')
    plt.title('Proportion of Target Group in Each Salinity Category')
    
    #save plot
    output_plot_path = str(output_file + "_salinity.png")
    plt.savefig(output_plot_path, bbox_inches='tight')
    
    print("Plot saved to file.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the salinity range of your group of interest.')
    parser.add_argument('-s', '--samples', type=str, help='Path to samples file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    main(args.samples, args.output)