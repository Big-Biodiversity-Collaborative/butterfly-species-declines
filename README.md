# Butterfly Species Declines

Code and infrastructure for building temporally-specific (pre- and post-2000) 
species distribution models for a variety of butterfly species. The project is 
linked to work and data in a long-term transect of butterfly observations 
across CA and NV, with the ultimate goal of providing infrastructure to make 
comparisons of distributions of all butterfly species in North America between
the two time periods.  

The project can be split into several components:  
1. [Automated extraction of occurrence records for species from iNat and GBIF](./R/pulling_organizing_data)  
2. Prepping data pre-modeling (combining environmental rasters, cleaning, 
background point generation, splitting data into training and testing data).  
3. Model building  
4. Model evaluation  
5. Predictions, maps and visualizations  

## General overview of major modeling steps  
The steps outlined below are automated in the `build_sdm()` function and 
outputs a master list of objects that are the complete step-wise record of the
following functions.  

1. The `prep_data()` function prepares individual occurence data for further 
analysis by splitting it into the two time periods (T1, T2), generating 
10k background points, cropping the environmental rasters to the occurence data
and generating a SpatialPointsDataFrame object that is used in further 
analysis.  

2. The `run_block_cv()` function makes use of the [blockCV](https://github.com/rvalavi/blockCV) package developed by R Valavi. 
Here, we use build a spatial grid across the area of interest and split the 
data into training (blocks 1-4) and test data (block 5) that we use to evaluate 
models.  

3. The `prep_data_2()` function does anothe round of data manipulation, this 
time extracting the environmental data for the occurence and background points, 
filtering out some unwanted columns and dropping missing data and filtering out 
some bad data (-Inf).  

4. The `train_test_split()` function takes occurence data and parses it into 
training and test sets, based on the methodology set up by `run_block_cv()`.  



There are also resources to run analyses in parallel UArizona's HPC (e.g. PBS 
scripts). 
