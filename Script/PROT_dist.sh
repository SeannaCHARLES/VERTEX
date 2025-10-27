#!/bin/bash

# This file launches DNAdist and create a distance matrix between the organismes. 
# Using the organismes ID as row and column header

##########################Input##########################
path="$1"
allstudies="$2"
resfile="$3"
IFS=' ' read -ra studies <<< "$allstudies"
##########################################################

###########################Script##########################
for study in "${studies[@]}"; do 

        echo A DNAdist table is created for the "${study}" study 
#/home/seannacharles/Peribunyaviridae/res_phylo/NSmGn_Peribunyaviridae_nucl_OG.aln-fasta
        seqret -sequence "${path}${study}${resfile}NSmGn_${study}_prot_OG.aln-fasta" -outseq "${path}${study}${resfile}*${study}_prot.phy" -osformat phylip 
        #converts a fasta file into phylip
        
        cp "${path}${study}${resfile}*${study}_prot.phy" infile
        #phylip DNAdist need a phylip format as input. The input need to be named infile

        echo -e "Y\n" | phylip protdist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read -r i; do 
        #All the id supposedly start with a letter this caracteristique is used to identify the organismes ID

        echo -ne "\t${i}"
        done # for the while loop
        echo "" ;} > "${path}${study}/${study}_protdist.txt"


        nb_colonnes=$(head -n1 "${path}${study}/${study}_protdist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max="$nb_colonnes" '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${study}/${study}_protdist.txt"
        rm outfile >/dev/null
        rm infile >/dev/null

awk 'NR==1 {                # Première ligne = en-tête
    printf "%s", $1          # Affiche la 1re colonne (nom des séquences)
    for (i=2; i<=NF; i++)    # Parcourt les colonnes
        if ($i ~ /^N/) keep[i]=1
        else keep[i]=0
    for (i=2; i<=NF; i++) if (keep[i]) printf "\t%s", $i
    printf "\n"
    next
}
$1 ~ /^G/ {                  # Garde les lignes dont le nom commence par G
    printf "%s", $1
    for (i=2; i<=NF; i++) if (keep[i]) printf "\t%s", $i
    printf "\n"
}' "${path}${study}/${study}_protdist.txt" > "${path}${study}/${study}_protdist_filter.txt"

done # for the for loop
############################################################