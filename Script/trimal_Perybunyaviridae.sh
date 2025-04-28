#!/bin/bash

path="/home/seannacharles/Documents/Perybunyaviridae/seg_M_phylo/"
pwd

Type_seq=("nucl" "prot")
for type in "${Type_seq[@]}"; do
        File=("Gc${type}" "Gn${type}" "NSm${type}" "NSmGn${type}")
        for ID in "${File[@]}"; do 
                echo "${ID} FASTA"
                trimal  -in "${path}${ID}_align.aln-fasta" \
                        -out "${path}${ID}_trimal.fasta" \
                        -automated1 \
                        -fasta 
                echo "${ID} PHYLIP"
                trimal  -in "${path}${ID}_align.aln-fasta" \
                        -out "${path}${ID}_trimal.phy" \
                        -automated1 \
                        -phylip 
        done
done



File=("Gcnucl" "Gnnucl" "NSmnucl" "NSmGnnucl")
for ID in "${File[@]}"; do 
        echo dnadist 
        cp "${path}${ID}_trimal.phy" infile
        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
        echo -ne "\t${i}"
        done
        echo "" ;} > "${path}${ID}_dnadist.txt"


        nb_colonnes=$(head -n1 "${path}${ID}_dnadist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max="$nb_colonnes" '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  # Nouvelle ligne après 10 valeurs
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${ID}_dnadist.txt"
        rm outfile
done


study="Perybunyaviridae" 
ths=0.01
target="bunyamwer"
proteins=("Gn" "NSm" "Gc")
list_proteins=$(IFS=, ; echo "${proteins[*]}")

/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths" "$list_proteins"
