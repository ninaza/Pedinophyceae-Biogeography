#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:59:55 2024

@author: ninpo556
"""

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from shapely.geometry import Point

# Load sample data (replace with your file path)
samples = pd.read_csv("data/raw/eukbank_18S_V4_samples.tsv", sep='\t')  # Columns: 'Sample', 'Longitude', 'Latitude'
sample_lon_lat = samples[["sample", "latitude", "longitude"]].copy()

# Convert sample data into GeoDataFrame
geometry = [Point(xy) for xy in zip(sample_lon_lat['longitude'], samples['latitude'])]
samples_gdf = gpd.GeoDataFrame(sample_lon_lat, geometry=geometry, crs="EPSG:4326")

# Load Longhurst province polygons (replace with your shapefile or GeoJSON file path)
longhurst = gpd.read_file("data/raw/GOaS_v1_20211214/goas_v01.shp")  # Replace with actual file

# Ensure CRS match
longhurst = longhurst.to_crs(samples_gdf.crs)

# Perform spatial join to assign Longhurst codes to samples
samples_with_longhurst = gpd.sjoin(samples_gdf, longhurst, how="left", predicate="within")

# Save the result to a CSV file
samples_with_longhurst.to_csv("data/edited/eukbank_18S_V4_samples_with_longhurst.csv", index=False)
