#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import argparse

def main(input_file, output_file):

    sample_df = pd.read_csv(input_file, sep='\t', usecols=["latitude", "longitude"])
    
# initialize an axis
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

# plot map on axis
    countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    countries.plot(color="lightgrey", ax=ax)

# plot points
    sample_df.plot(x="longitude", y="latitude", kind="scatter", ax=ax)

    plt.savefig(output_file, bbox_inches='tight')

    print("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process TSV file and generate PNG file.')
    parser.add_argument('input_file', type=str, help='Path to input TSV file')
    parser.add_argument('output_file', type=str, help='Path to output PNG file')
    args = parser.parse_args()

    main(args.input_file, args.output_file)
