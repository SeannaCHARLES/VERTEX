#!/bin/bash


##CCHFV
Type_seq=("nucl")
path="/home/seannacharles/Documents/CCHFV/seg_M_phylo/"
pwd



for type in "${Type_seq[@]}"; do
        echo "conversion des sequences ${type}"
        File=("Gc${type}" "Gn${type}" "GP38${type}" "MLD${type}" "NSm${type}" "NSmGn${type}")
        for ID in "${File[@]}"; do 
                seqret -sequence "${path}${ID}.aln-fasta" -outseq "${path}${ID}_trimal.phy" -osformat phylip
                echo dnadist 
                cp "${path}${ID}_trimal.phy" infile
                echo -e "Y\n" | phylip dnadist
                { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
                echo -ne "\t${i}"
                done
                echo "" ;} > "${path}${ID}_dnadist.txt"


                nb_colonnes=$(head -n1 "${path}${ID}_dnadist.txt" | awk -F"\t" '{print NF}')

                sed '1d' outfile |awk -v max=$nb_colonnes '{
                for (i=1; i<=NF; i++) {
                printf "%s\t", $i;
                count++;
                if (count % max == 0) print "";  # Nouvelle ligne après 10 valeurs
                }
                }nucl
                END { if (count % max != 0) print ""; }' >> "${path}${ID}_dnadist.txt"
                rm outfile
        done

done



echo Prottest CCHFV


study="CCHFV"     
target="np950235."
ths=0.009

/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths"
