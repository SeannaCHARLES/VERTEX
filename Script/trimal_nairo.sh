#!/bin/bash


echo Trimal Nairoviridae 
path="/home/seannacharles/Documents/Nairoviridae/seg_M_phylo/"
pwd
Type_seq=("nucl" "prot")
for type in "${Type_seq[@]}"; do
        echo "Traitement des sequence ${type}"
                File=("Gc${type}" "Gn${type}" "NSm${type}" "NSmGn${type}" "GP38${type}" "MLD${type}")
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

File=("Gcnucl" "Gnnucl" "NSmnucl" "NSmGnnucl" "MLDnucl" "GP38nucl")
for ID in "${File[@]}"; do 
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
        }
        END { if (count % max != 0) print ""; }' >> "${path}${ID}_dnadist.txt"
        rm outfile
done


study="Nairoviridae"     
target="haemorrha"
ths=0.005
# Convertir la liste en une chaîne séparée par des virgules


/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths"
# echo prottest nairo


 

# echo NSmGn : >/home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/NSmGn_trimal.fasta -all -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# echo NSm : >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/NSm_trimal.fasta -all -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# echo Gn : >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Gn_trimal.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# echo Gc : >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Gc_trimal.fasta -all  -DT -F  |grep "Best model according to DT:">> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# echo GP38 : >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/GP38_trimal.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# echo MLD : >> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt
# prottest -i /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/MLD_trimal.fasta -all  -DT -F  |grep "Best model according to DT:">> /home/seannacharles/Documents/Nairoviridae/seg_M_phylo/Best_model_nairo.txt

# sudo  mv /usr/bin/phyml /usr/bin/phyml_3.1 
# sudo  mv /usr/bin/phyml_3.3 /usr/bin/phyml 