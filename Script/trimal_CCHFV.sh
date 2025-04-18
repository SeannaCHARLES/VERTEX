#!/bin/bash


##CCHFV
Type_seq=("nucl")
echo Trimal CCHFV 
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
# Convertir la liste en une chaîne séparée par des virgules


/usr/bin/python3 /home/seannacharles/Documents/Script/delete_dbl.py "$study"  "$target" "$ths"
# sudo  mv /usr/bin/phyml /usr/bin/phyml_3.3 
# sudo  mv /usr/bin/phyml_3.1 /usr/bin/phyml 
 
# echo NSmGn : > /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/NSmGn_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# echo NSm : >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt
# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/NSmpro_align_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# echo MLD : >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt
# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/MLD_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# echo GP38 : >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt
# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/GP38_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# echo Gn : >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt
# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/Gn_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# echo Gc : >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt
# prottest -i /home/seannacharles/Documents/CCHFV/seg_M_phylo/Gc_trimal_sdbl.fasta -all  -DT -F |grep "Best model according to DT:" >> /home/seannacharles/Documents/CCHFV/seg_M_phylo/Best_model_CCHFV.txt

# sudo  mv /usr/bin/phyml /usr/bin/phyml_3.1 
# sudo  mv /usr/bin/phyml_3.3 /usr/bin/phyml 