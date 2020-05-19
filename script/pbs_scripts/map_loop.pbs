#!/bin/bash

# NOTE: this script is an example - it worked on a smaller version of this project
# when things were structured different. It will need to be reconfigured for the
# newer dev version.

#PBS -N map_loop
#PBS -W group_list=klprudic
#PBS -q standard
#PBS -l select=1:ncpus=16:mem=96gb:pcmem=6gb
#PBS -l walltime=24:0:0
#PBS -j oe

# Loading modules on the HPC
module load java/1.8 proj/5/5.2.0 gdal/2/2.4.2 geos/3.5/3.5.1 R

# Getting rJava to work
R CMD javareconf -e
export R_JAVA_LD_LIBRARY_PATH=${JAVA_HOME}/jre/lib/amd64/server

# Running the mapping loop

cd /extra/keatonwilson/butterfly-species-declines/
Rscript ./script/analysis_mapping/map_loop.R
