#!/bin/bash

path="/home/seannacharles/Documents/"
study="Nairoviridae"
phylofile="/seg_M_phylo/"

Type_seq=("nucl" "prot")
File=("Gc${type}" "Gn${type}" "NSm${type}" "NSmGn${type}" "GP38${type}" "MLD${type}")

for type in "${Type_seq[@]}"; do
        echo "Traitement des sequence ${type}"
        
        for ID in "${File[@]}"; do 
                echo "${ID} FASTA"
                trimal  -in "${path}${study}${phylofile}${ID}_align.aln-fasta" \
                        -out "${path}${study}${phylofile}${ID}_trimal.fasta" \
                        -automated1 \
                        -fasta 
                echo "${ID} PHYLIP"
                trimal  -in "${path}${study}${phylofile}${ID}_align.aln-fasta" \
                        -out "${path}${study}${phylofile}${ID}_trimal.phy" \
                        -automated1 \
                        -phylip 
                        
        done
done

File=("Gcnucl" "Gnnucl" "NSmnucl" "NSmGnnucl" "MLDnucl" "GP38nucl")

for ID in "${File[@]}"; do 
        echo dnadist 
        cp "${path}${study}${phylofile}${ID}_trimal.phy" infile
        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
        echo -ne "\t${i}"
        done
        echo "" ;} > "${path}${study}${phylofile}${ID}_dnadist.txt"


        nb_colonnes=$(head -n1 "${path}${study}${phylofile}${ID}_dnadist.txt" | awk -F"\t" '{print NF}')

        sed '1d' outfile |awk -v max=$nb_colonnes '{
        for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  # Nouvelle ligne après 10 valeurs
        }
        }
        END { if (count % max != 0) print ""; }' >> "${path}${study}${phylofile}${ID}_dnadist.txt"
        rm outfile >/dev/null
done


study="Nairoviridae"     
target="haemorrha"
ths=0.005
proteins=("MLD" "GP38" "Gn" "NSm" "Gc")
list_proteins=$(IFS=, ; echo "${proteins[*]}")


/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths" "$list_proteins"
