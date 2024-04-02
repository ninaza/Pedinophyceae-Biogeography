#!/usr/bin/env python 

####### DISTRIBUTION AND ABUNDANCE ############################################

import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs

### PLOTTING SAMPLING POINTS ###
# open group sample table
sample_df = pd.read_csv("/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/parserTest/Pedino_samples.tsv", sep='\t')

# extract longitude and latitude data
lon_lat_df = sample_df[["latitude", "longitude"]].copy()

# initialize an axis
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

# plot map on axis
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
countries.plot(color="lightgrey", ax=ax)

# plot points
lon_lat_df.plot(x="longitude", y="latitude", kind="scatter", ax=ax)

# save plot to file
plt.savefig('/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.samples.png', bbox_inches='tight')

### PLOTTING SAMPLING POINTS INCLUDING ABUNDANCE INFORMATION ###
# open counts table
counts_df = pd.read_csv("/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.counts.tsv", sep='\t')

## calculation of relative abundance of Pedinophyceae endosymbionts per sample
# Group by 'Sample' and calculate the mean number of reads per sample
average_reads_per_sample = counts_df.groupby('sample')['nreads'].mean()
average_reads_per_sample = average_reads_per_sample.rename('nreads-Pedino')

# take sample and nreads columns from sample_df, set sample name as index and rename nreads to nreads total
total_sample_abundance_df = sample_df[["sample", "nreads"]].copy()
total_sample_abundance_df = total_sample_abundance_df.set_index('sample')
total_sample_abundance_df = total_sample_abundance_df.rename(columns={'nreads': 'nreads-total'})

# merge df with mean abundance of group per sample and df with total number of reads in sample and calulcate ratio and percentage
merged_df = total_sample_abundance_df.merge(average_reads_per_sample, left_on='sample', right_index=True)
merged_df['relative_abundance'] = merged_df['nreads-Pedino'] / merged_df['nreads-total']
merged_df['abundance_percent'] = merged_df['relative_abundance'] * 100

# save merged df to file
output_file_path = "/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/parserTest/Pedino-EukBank_nreads.tsv"
merged_df.to_csv(output_file_path, sep='\t')

# get max abundance value and name of sample
max_abundance = merged_df.loc[merged_df['relative_abundance'].idxmax()]
max_abundance_sample = max_abundance.name
max_abundance_sample_info = sample_df[sample_df['sample'] == max_abundance_sample]
print(max_abundance_sample)

## make a map plot with dot size representing abundance
lons = sample_df['longitude'].tolist()
lats = sample_df['latitude'].tolist()
percentages = merged_df['abundance_percent'].tolist()

# Define the maximal value for abundance
#max_abundance = max(percentages)
max_abundance = 10

# Calculate the percentage values based on the maximal value
percentage_values = [max_abundance, max_abundance / 2, max_abundance / 4, max_abundance / 8]

# Create a new figure
plt.figure(figsize=(12, 8), dpi=300)

# Create a GeoAxes with the PlateCarree projection
ax = plt.axes(projection=ccrs.PlateCarree())

# Add coastlines
ax.coastlines()

# Plot the data with varying marker sizes based on percentage values
scatter = ax.scatter(lons, lats, s=[p * 100 for p in percentages], c='green', alpha=0.5, transform=ccrs.PlateCarree())

# Create legend scatter plots for reference sizes
legend_handles = []
for pct in percentage_values:
    size = pct * 100  # Scaling the size for better visualization
    handle = ax.scatter([], [], s=size, label=f'{pct:.2f}%', c='green', alpha=0.5)
    legend_handles.append(handle)

# Add legend
ax.legend(handles=legend_handles, title='Relative abundance', loc='lower left', bbox_to_anchor=(1, 0.5),)

#save plot to file
plt.savefig('/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.abundance.png',
            bbox_inches = 'tight')

######### SIZE FRACTION #######################################################

max_lower_threshold = sample_df.loc[sample_df['size_fraction_lower_threshold'].idxmax()]
samples_max_lower = sample_df[sample_df['size_fraction_lower_threshold'] == 180]
min_upper_threshold = sample_df.loc[sample_df['size_fraction_upper_threshold'].idxmin()]

sample_names = sample_df['sample'].tolist()
size_fraction_lower = sample_df['size_fraction_lower_threshold'].tolist()
size_fraction_upper = sample_df['size_fraction_upper_threshold'].tolist()

# Create a DataFrame to store sample names and depths
size_fraction_df = pd.DataFrame({'Sample': sample_names, 'Size fraction lower µm': size_fraction_lower, 'Size fraction upper µm': size_fraction_upper})

size_bins = [-float('inf'), 1.9, 20, 200, float('inf')]
size_labels = ['pico', 'nano', 'micro', 'meso']
size_fraction_df['Size category lower'] = pd.cut(size_fraction_df['Size fraction lower µm'], bins=size_bins, labels=size_labels)
size_fraction_df['Size category upper'] = pd.cut(size_fraction_df['Size fraction upper µm'], bins=size_bins, labels=size_labels)


######## DEPTH ################################################################

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

#save plot to file
plt.savefig('/Users/ninpo556/Mirror/ASPP-course/Pedinophyceae-Biogeography/Pedino-EukBank.depth.png',bbox_inches = 'tight')


######## SALINITY #############################################################

######## BIOME ################################################################