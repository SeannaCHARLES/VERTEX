#!/bin/bash


path="/home/seannacharles/Documents/"
study="CCHFV"
phylofile="/seg_M_phylo/"

File=("Gcnucl" "Gnnucl" "GP38nucl" "MLDnucl" "NSmnucl" "NSmGnnucl")

#conversion des sequences nucl

for ID in "${File[@]}"; do 
        seqret -sequence "${path}${study}${phylofile}${ID}.aln-fasta" -outseq "${path}${study}${phylofile}${ID}_trimal.phy" -osformat phylip
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
                        if (count % max == 0) print "";  
                        }
                        }nucl
                        END { if (count % max != 0) print ""; }' >> "${path}${study}${phylofile}${ID}_dnadist.txt"
        rm outfile >/dev/null
done




   
target="np950235."
ths=0.009
proteins=("MLD" "GP38" "Gn" "NSm" "Gc")
list_proteins=$(IFS=, ; echo "${proteins[*]}")

/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths" "$list_proteins"
