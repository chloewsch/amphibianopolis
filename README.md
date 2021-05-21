# amphibianopolis

Code accompanying Schmidt & Garroway (2021) The population genetics of urban and rural amphibians in North America

Data:
amph_data.csv: final dataset
- pop: unique site identifier	
- species: formatted as Genus_species	
- author: author of original data source (see Table S2 for full list of references)	
- order: taxonomic order; Anura/Caudata	
- lon, lat: longitude & latitude coordinates in decimal degrees, WGS84
- num_individuals: number of individuals sampled
- gene_diversity: Nei's gene diversity
- allelic_richness: rarefied allelic richness (min. 5 individuals)	
- global_fst: population-specific FST
- Ne: effective population size estimate
- Ne_upper,	Ne_lower: upper and lower bounds of effective population size	
- num_loci: number of loci sampled
- urban: 1 = urban sample site; 0 = non-urban sample site	
- road_1_km; road_5_km;	road_10_km;	road_15_km: road lengths within 1, 5, 10, 15 km buffers around sites
- pd_1_km; pd_5_km; pd_10_km; pd_15_km: mean human population density within 1, 5, 10, 15 km buffers around sites 
- hfi_1_km;	hfi_5_km; hfi_10_km; hfi_15_km: mean Human Footprint Index within 1, 5, 10, 15 km buffers around sites
- imp_30m; imp_1_km; imp_5_km; imp_10_km; imp_15_km: mean % impervious surface cover within 30m and 1, 5, 10, 15 km buffers around sites


R code:
1. build_urban_data.R
- scripts to extract urban variables & create final dataset used in analyses

2. 2_dbmem_regressions_scales.R
- compute, then select significant MEMs
- Bayesian GLMMs for all variables across scales
- plots

3. 3_phylo_regressions
- models and plots for phylogenetic regressions (main text models only)
