#!/bin/bash

# This file launches DNAdist and create a distance matrix between the organismes. 
# Using the organismes ID as row and column header

##########################Input##########################
path="$1"
allstudies="$2"
IFS=' ' read -ra studies <<< "$allstudies"
##########################################################

###########################Script##########################
for study in "${studies[@]}"; do 

        echo A DNAdist table is created for the "${study}" study 

        seqret -sequence "${path}${study}/nucl_all_align_${study}.aln-fasta" -outseq "${path}${study}/nucl_all_align_${study}.phy" -osformat phylip 
        #converts a fasta file into phylip
        
        cp "${path}${study}/nucl_all_align_${study}.phy" infile
        #phylip DNAdist need a phylip format as input. The input need to be named infile

        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read -r i; do 
        #All the id supposedly start with a letter this caracteristique is used to identify the organismes ID

        echo -ne "\t${i}"
        done # for the while loop
        echo "" ;} > "${path}${study}/${study}_dnadist.txt"


        nb_colonnes=$(head -n1 "${path}${study}/${study}_dnadist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max="$nb_colonnes" '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${study}/${study}_dnadist.txt"
        rm outfile >/dev/null
        rm infile >/dev/null
done # for the for loop
############################################################