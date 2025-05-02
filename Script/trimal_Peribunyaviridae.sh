#!/bin/bash

study="$1"
Type_seq=("prot" "nucl")
for type in "${Type_seq[@]}"; do
        File=("Gc${type}" "Gn${type}" "NSm${type}" "NSmGn${type}")
        for ID in "${File[@]}"; do 
                echo "${ID} FASTA"
                trimal  -in "${path}${study}${res_file}${ID}_align.aln-fasta" \
                        -out "${path}${study}${res_file}${ID}_trimal.fasta" \
                        -automated1 \
                        -fasta \
                        -block 5
                echo "${ID} PHYLIP"
                trimal  -in "${path}${study}${res_file}${ID}_align.aln-fasta" \
                        -out "${path}${study}${res_file}${ID}_trimal.phy" \
                        -automated1 \
                        -phylip \
                        -block 5
        done
done



File=("Gcnucl" "Gnnucl" "NSmnucl" "NSmGnnucl")

for ID in "${File[@]}"; do 
        echo dnadist 
        cp "${path}${study}${res_file}${ID}_trimal.phy" infile
        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
        echo -ne "\t${i}"
        done
        echo "" ;} > "${path}${study}${res_file}${ID}_dnadist.txt"


        nb_colonnes=$(head -n1 "${path}${study}${res_file}${ID}_dnadist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max="$nb_colonnes" '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${study}${res_file}${ID}_dnadist.txt"
        rm outfile 
        rm infile
done



ths=0.01
target="bunyamwer"
proteins=("Gn" "NSm" "Gc")
list_proteins=$(IFS=, ; echo "${proteins[*]}")

/usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "$study"  "$target" "$ths" "$list_proteins"
