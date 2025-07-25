#!/bin/bash
SECONDS=0

while IFS=$',' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 _; do
    # Ignorer les lignes commentées
    [[ $col1 =~ ^# ]] && continue
    names+=("$col2") 
    Allstudies+=("$col1")
    for_prot+=("$col5")

    if [[ $col4 == N* ]]; then
        studies+=("$col1")  #Studies associé aux names 
        NUID+=("$col3") #Identifiants des outgroup spécifique
        target_list+=("$col8")  #Organisme représentant l'étude
        position+=("$col9") #Positions des proteines de l'outgroup
        start+=("$col10")   #Début de la séquence proteiques sur la séquence nucleique
        ths+=("$col11") #Seuil pour DNAdist pour supprimer les doublons
        proteins+=("$col6")
        if [[ $col14 == Y* ]]; then
            trim+=("$col1")
        fi
    else
        studiesspeOG+=("$col1") #Studies associé aux names 
           #Nom des 'sous unités' de chaque studies ex. spécifique outgroup
        NUIDspeOG+=("$col3")    #Identifiants des outgroup spécifique
        OG_prot+=("$col5")   #Proteines associées aux outgroup specifique
        positionspeOG+=("$col9")    #Positions des proteines de l'outgroup
        startspe_OG+=("$col10")  #Début de la séquence proteiques sur la séquence nucleique
        proteinsspe_OG+=("$col6")
    fi
done < All_information.txt

#echo "${NUID[*]}"
#NUID=$(grep -v "^#" All_information.txt |awk -F"," '{print $3}')
path=$(grep '^#' All_information.txt | awk -F"=" 'NR==2 {print $2}'| tr -d ',')
res_file=$(grep '^#' All_information.txt | awk -F"=" 'NR==3 {print $2}'| tr -d ',')
path_script=$(grep '^#' All_information.txt | awk -F"=" 'NR==4 {print $2}'| tr -d ',')
Type_seq=$(grep '^#' All_information.txt | awk -F"=" 'NR==5 {print $2}'| tr -d ',')
email=$(grep '^#' All_information.txt | awk -F"=" 'NR==6 {print $2}'| tr -d ',')
keep_prot=$(grep -v "^#" All_information.txt |awk -F"," '{print $12}')
#proteins=$(grep -v "^#" All_information.txt |awk -F"," '{print $6}'| grep -v '^["[:space:]"]*$'|sed 's/\,/" "/g'| uniq)
speprot=$(grep -v "^#" All_information.txt |awk -F"," '{print "\"",$7,"\""}')
orgaWIWO=$(grep -v "^#" All_information.txt | awk -F"," '{gsub(/\xc2\xa0/, "", $13); print $13}' | awk -F"=" '{print $2}')
INOUT=$(grep -v "^#" All_information.txt |awk -F"," '{print $13}'|awk -F"=" '{print $1}'| awk -F "*" '{print $1}')
protINOUT=$(grep -v "^#" All_information.txt |awk -F"," '{print $13}'|awk -F"=" '{print $1}'| awk -F "*" '{print $2}')



for study in "${studies[@]}"; do
rm -r "${path}""${study:?}/"* 
mkdir -p "${path}""${study}"/nucl "${path}""${study}"/prot "${path}""${study}""${res_file}"/treefolder "${path}""${study}"/legend "${path}""${study}"/blast/database/
done

# #Common scripts for all studies
 echo The sequences are currently being downloaded ...



for i in "${!NUID[@]}"; do   
    study="${studies[$i]}"
    echo "${study}"
    IFS=$';' read -ra IDprot <<< "${NUID[$i]}"
   for ID in "${IDprot[@]}"; do 

        esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb
        cat "${path}""${study}"/"$ID".gb >> "${path}""${study}"/sequences_"${study}".gb
        rm "${path}""${study}"/"$ID".gb >/dev/null
    done
done


for i in "${!NUIDspeOG[@]}"; do   
    study="${studiesspeOG[$i]}"
    IFS=$';' read -ra IDprotspeOG <<< "${NUIDspeOG[$i]}"
    for ID in "${IDprotspeOG[@]}"; do 

        esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb

    done
done



echo The captions, metadata files and files containing all the sequences are currently being created.
/usr/bin/python3 "${path}""${path_script}"genbank_extract.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${names[*]}" "${NUIDspeOG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${email}"

./"${path_script}"DNA_dist.sh "${path}" "${studies[*]}" 

echo The duplicate sequences are being deleted
/usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "${studies[*]}" "${keep_prot[*]}" "${ths[*]}" "${target_list[*]}" 


echo the proteins are being seperated
/usr/bin/python3 "${path}""${path_script}"test_script_mark_spe.py "${path}" "${studies[*]}" "${position[*]}" "${names[*]}" "${positionspeOG[*]}" "${start[*]}" "${startspe_OG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${res_file}" "${studiesspeOG[*]}" "${NUIDspeOG[*]}"

###Blast pour avoir le % de similarité table pour chaque study
./"${path_script}"Blast_prot.sh "${path}" "${res_file}" "${studies[*]}" 
for study in "${studies[@]}"; do
    echo "${study}" blast
        for fichier in "${path}${study}${res_file}"*prot_OG.fasta; do
            nom_fichier=$(basename "$fichier")
            start_file_name="${nom_fichier%prot_OG.fasta}"

            makeblastdb -in "${path}${study}${res_file}${start_file_name}_OG.fasta"\
                        -title "${start_file_name}"\
                        -dbtype prot \
                        -out "${path}${study}/blast/database/${start_file_name}_DB"\
                        -parse_seqids

        done

        for fichier in "${path}${study}"/blast/database/*_DB; do
            nom_fichier=$(basename "$fichier")
            start_file_name="${nom_fichier%_DB}"

            blastp  -query "${path}${study}${res_file}${start_file_name}_OG.fasta" \
                    -db "${fichier}" \
                    -out "${path}${study}/blast/${start_file_name}_blast.txt" \
                    -outfmt 7 

        done
done

echo The proteines that are presented or not in certains proteins are being treated
/usr/bin/python3 "${path}""${path_script}"WIWO_protein.py "${path}" "${names[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${res_file}" "${Type_seq[*]}"

#merge = OG

/usr/bin/python3 "${path}""${path_script}"merge_spe_prot.py "${path}" "${names[*]}" "${speprot[*]}" "${Type_seq[*]}" "${res_file}" "${studies[*]}" "${Allstudies[*]}" "${for_prot[*]}" 


for study in "${studies[@]}"; do
    echo "${study}"
    for val in "${trim[@]}"; do
        for fichier in "${path}${study}${res_file}"*_OG.aln-fasta; do
            nom_fichier=$(basename "$fichier")
            start_file_name="${nom_fichier%_OG.aln-fasta}"
            if [[ "$study" == "$val" ]]; then


                    trimal  -in "$fichier" \
                        -out "${path}${study}${res_file}${start_file_name}/treefolder/_OG_trimal.fasta" \
                        -automated1 \
                        -fasta 
                    trimal  -in "$fichier" \
                        -out "${path}${study}${res_file}${start_file_name}/treefolder/__OG_trimal.phy" \
                        -automated1 \
                        -phylip 
            else 
                cp "$fichier" "${path}${study}${res_file}${start_file_name}/treefolder/__OG_trimal.fasta"
                seqret -sequence "$fichier" -outseq "${path}${study}${res_file}${start_file_name}/treefolder/__OG_trimal.phy" -osformat2 phylip
            { echo "${start_file_name}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${study}""${res_file}""${start_file_name}"_OG_trimal.phy ;}>> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            iqtree2 -s "${fichier}" -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            echo "\n" >>"${path}""${study}""${res_file}"Best_model_"${study}".txt
            fi
        done
    done
done


# for study in "${studies[@]}"; do
#     for fichier in "${path}${study}${res_file}"*_OG_trimal.fasta; do
#         nom_fichier=$(basename "$fichier")
#         start_file_name="${nom_fichier%_OG_trimal.fasta}"    
#         { echo "${start_file_name}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${study}""${res_file}""${start_file_name}"_OG_trimal.phy ;}>> "${path}""${study}""${res_file}"Best_model_"${study}".txt
#         iqtree2 -s "${fichier}" -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}""${res_file}"Best_model_"${study}".txt

#     done

# done    


Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "It is finally over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"