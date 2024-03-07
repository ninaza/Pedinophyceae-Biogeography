#!/usr/bin/env python 

import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd

group_sample_df = pd.read_csv("/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.tsv", 
                              sep='\t', usecols=["latitude", "longitude"])
# initialize an axis
fig, ax = plt.subplots(figsize=(8,6))

# plot map on axis
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
countries.plot(color="lightgrey", ax=ax)

# plot points
group_sample_df.plot(x="longitude", y="latitude", kind="scatter", ax=ax)

plt.savefig('/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.png', 
            bbox_inches='tight')
