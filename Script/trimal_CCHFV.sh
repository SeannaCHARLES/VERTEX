#!/bin/bash

study="$1"
path="$2"
res_file="$3"
path_script="$4"

# path="/home/seannacharles/"
# path_script="Phylogenetic/Script/"
# res_file="/res_phylo/"
# study="CCHFV"

File=("Gcnucl" "Gnnucl" "GP38nucl" "MLDnucl" "NSmnucl" "NSmGnnucl")


for ID in "${File[@]}"; do 
echo "${path}${study}${res_file}${ID}_align.aln-fasta"
        seqret -sequence "${path}${study}${res_file}${ID}_align.aln-fasta" -outseq "${path}${study}${res_file}${ID}_trimal.phy" -osformat phylip
        cp "${path}${study}${res_file}${ID}_trimal.phy" infile
        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
        echo -ne "\t${i}"
        done
        echo "" ;} > "${path}${study}${res_file}${ID}_dnadist.txt"


        nb_colonnes=$(head -n1 "${path}${study}${res_file}${ID}_dnadist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max=$nb_colonnes '{
                for (i=1; i<=NF; i++) {
                        printf "%s\t", $i;
                        count++;
                        if (count % max == 0) print "";  
                        }
                        }nucl
                        END { if (count % max != 0) print ""; }' >> "${path}${study}${res_file}${ID}_dnadist.txt"
        rm outfile 
        rm infile
done




   
target="np950235."
ths=0.009
proteins=("MLD" "GP38" "Gn" "NSm" "Gc")
list_proteins=$(IFS=, ; echo "${proteins[*]}")

/usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "$study"  "$target" "$ths" "$list_proteins" "$res_file"
