#!/bin/bash

#setting path to pbs starter
PBS_starter=../../script/pbs_scripts/pbs.tmp

# listing files in split data folder
cd ./data/split_data/
ls

for i in `ls ./`
do
        # print to standard_out the content of PBS_starter, find and replace 
        # every instance of "SPECIES" with the value of $i, this is written 
        # to a file (line by line)
        # the file will be overwritten if present 
        
        cat $PBS_starter | sed "s/SPECIES/$i/g" > ${PBS_starter%.*}.$i.pbs 
        #run that shit
        cd ../../output/pbs_out
        qsub ${PBS_starter%.*}.$i.pbs  
done	