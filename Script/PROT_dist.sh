#!/bin/bash

# This file launches DNAdist and create a distance matrix between the organismes. 
# Using the organismes ID as row and column header

##########################Input##########################
file="$1"
study="$2"
path="$3"

##########################################################

###########################Script##########################
#mafft --auto "${path}${file}.fasta > ${path}${file}.aln-fasta"
seqret -sequence "${path}${file}.aln-fasta" -outseq "${path}${file}.phy" -osformat phylip 
        #converts a fasta file into phylip     
cp "${path}${file}.phy" infile
        #phylip DNAdist need a phylip format as input. The input need to be named infile
echo -e "Y\n" | phylip protdist
{ sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read -r i; do 
        #All the id supposedly start with a letter this caracteristique is used to identify the organismes ID
echo -ne "\t${i}"
        done # for the while loop
echo "" ;} > "${path}${file}_protdist.txt"


nb_colonnes=$(head -n1 "${path}${file}_protdist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max="$nb_colonnes" '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${file}_protdist_.txt"
        rm outfile >/dev/null
        rm infile >/dev/null

# awk 'NR==1 {                # Première ligne = en-tête
#     printf "%s", $1          # Affiche la 1re colonne (nom des séquences)
#     for (i=2; i<=NF; i++)    # Parcourt les colonnes
#         if ($i ~ /^N/) keep[i]=1
#         else keep[i]=0
#     for (i=2; i<=NF; i++) if (keep[i]) printf "\t%s", $i
#     printf "\n"
#     next
# }
# $1 ~ /^G/ {                  # Garde les lignes dont le nom commence par G
#     printf "%s", $1
#     for (i=2; i<=NF; i++) if (keep[i]) printf "\t%s", $i
#     printf "\n"
# }' "${path}${file}_protdist.txt" > "${path}${file}_protdist_.txt"


# tail -n +5 "${path}${study}/legend/Label_${study}_itol.txt" | grep -v '^#' > "${path}${study}/map.txt"
# awk -v OFS="\t" 'NR==FNR { map[$1]=$2; next }       # lit le fichier de correspondance (map.txt)
# {
#   # pour chaque champ (colonne) de la ligne courante
#   for (i=1; i<=NF; i++) {
#     if ($i in map) $i = map[$i]         # remplace si correspondance trouvée
#   }
#   print                                 # imprime la ligne modifiée
# }' "${path}${study}/map.txt" "${path}${file}_protdist_.txt"  >> "${path}${file}_protdist_filter.txt"

 # for the for loop
############################################################