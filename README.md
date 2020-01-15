# Butterfly Species Declines

Code and infrastructure for building temporally-specific species distribution 
models for a variety of butterfly species. The project is linked to work and 
data in a long-term transect of butterfly observations across CA and NV.  

The project can be split into several components:  
1. [Automated extraction of occurrence records for species from iNat and GBIF](./R/pulling_organizing_data)  
2. Prepping data pre-modeling (combining environmental rasters, cleaning, 
background point generation, splitting data into training and testing data).  
3. Model building  
4. Model evaluation  
5. Predictions, maps and visualizations  

There are also resources to run analyses in parallel UArizona's HPC (e.g. PBS 
scripts). 
