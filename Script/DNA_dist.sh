#!/bin/bash
path="$1"
allstudies="$2"
echo "${studies[@]}"
IFS=' ' read -ra studies <<< "$allstudies"

for study in "${studies[@]}"; do 
        echo dnadist 
        seqret -sequence "${path}${study}/nucl_all_align_${study}.aln-fasta" -outseq "${path}${study}/nucl_all_align_${study}.phy" -osformat phylip
        cp "${path}${study}/nucl_all_align_${study}.phy" infile
        echo -e "Y\n" | phylip dnadist
        { sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read -r i; do
        echo -ne "\t${i}"
        done
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
        rm outfile 
        rm infile
done