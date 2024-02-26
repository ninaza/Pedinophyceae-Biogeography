# Pedinophyceae-Biogeography
A biogeography study of Pedinophyceae using the EukBank metabarcoding database

In the coding project I will start working on extracting metadata (e.g. sampling site, sampling depth (ocean), sampling fraction) for my sequences of interest that belong to the a family of algae called Pedinophyceae. Extracting the data will not be very straight forward as each sequence can stem from different sampling efforst, meaning there are several metadatsets connected to each sequence. The metadat are also partially stored in different files. First I want to start collecting the data in a way that I can do data analysis on them later on. The data analysis will include some plotting (of geographical regions) as well as statistics.

- Counts table contains ASV names and sample names
- From there can make link to sample table and if necessary to project
- Find way to automatically check the project file for key-words like marine, arctic etc.

### Prepare sequences

- blasting again against EukBank only with the Pedinophyceae Blanes ASVs
- putting sequences in the Chlorophyte tree after to reassure it's only Pedinophyceae of the endosymbiont clade
- extract headers of the sequences
- start collecting the metadata
