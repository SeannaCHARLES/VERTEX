#!/bin/bash
SECONDS=0
########################################EXTRACT ALL THE INFORMATION FROM THE INFORMATION FILE############################################
awk -F',' '$0 !~ /^#/ && $13 == "" { $13 = "none*nn=nn" } { print }' OFS=',' All_information.txt > tmp && mv tmp All_information.txt
while IFS=$',' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 _; do
    
    [[ $col1 =~ ^# ]] && continue
    names+=("$col2") 
    Allstudies+=("$col1")
    for_prot+=("$col5")
    allID+=("$col3")
    
    if [[ $col4 == N* ]]; then
        studies+=("$col1")  #Studies associé aux names 
        NUID+=("$col3") #Identifiants des outgroup spécifique
        target_list+=("$col8")  #Organisme représentant l'étude
        position+=("$col9") #Positions des proteines de l'outgroup
        start+=("$col10")   #Début de la séquence proteiques sur la séquence nucleique
        ths+=("$col11") #Seuil pour DNAdist pour supprimer les doublons
        proteins+=("$col6")
        namespepro+=("$col7")
        if [[ $col14 == Y* ]]; then
            trim+=("$col1")
        fi
        diff_protname+=("$col15")
    else
        studiesspeOG+=("$col1") #Studies associé aux names 
           #Nom des 'sous unités' de chaque studies ex. spécifique outgroup
        NUIDspeOG+=("$col3")    #Identifiants des outgroup spécifique
        OG_prot+=("$col5")   #Proteines associées aux outgroup specifique
        positionspeOG+=("$col9")    #Positions des proteines de l'outgroup
        startspe_OG+=("$col10")  #Début de la séquence proteiques sur la séquence nucleique
        proteinsspe_OG+=("$col6")
        diff_protname_OG+=("$col15")
        namespeproOG+=("$col7")
    fi
done < All_information.txt

#echo "${NUID[*]}"
#NUID=$(grep -v "^#" All_information.txt |awk -F"," '{print $3}')
path=$(grep '^#' All_information.txt | awk -F"=" 'NR==2 {print $2}'| tr -d ',')
res_file=$(grep '^#' All_information.txt | awk -F"=" 'NR==3 {print $2}'| tr -d ',')
path_script=$(grep '^#' All_information.txt | awk -F"=" 'NR==4 {print $2}'| tr -d ',')
Type_seq=$(grep '^#' All_information.txt | awk -F"=" 'NR==5 {print $2}'| tr -d ',')
email=$(grep '^#' All_information.txt | awk -F"=" 'NR==6 {print $2}'| tr -d ',')
order=$(grep '^#' All_information.txt | awk -F":" 'NR==7 {print $2}'| tr -d ',')
keep_prot=$(grep -v "^#" All_information.txt |awk -F"," '{print $12}')
#proteins=$(grep -v "^#" All_information.txt |awk -F"," '{print $6}'| grep -v '^["[:space:]"]*$'|sed 's/\,/" "/g'| uniq)
speprot=$(grep -v "^#" All_information.txt |awk -F"," '{print "\"",$7,"\""}')
orgaWIWO=$(grep -v "^#" All_information.txt | awk -F"," '{gsub(/\xc2\xa0/, "", $13); print $13}' | awk -F"=" '{print $2}')
INOUT=$(grep -v "^#" All_information.txt |awk -F"," '{print $13}'|awk -F"=" '{print $1}'| awk -F "*" '{print $1}')
protINOUT=$(grep -v "^#" All_information.txt |awk -F"," '{print $13}'|awk -F"=" '{print $1}'| awk -F "*" '{print $2}')
diff_name_prot=$(grep "^#" All_information.txt |awk -F":" 'NR==7 {print $2}'| tr -d ',')
if [ ${#speprot[@]} -eq 1 ] && [ "${speprot[0]}" = '"  "' ]; then
speprot=()
fi

######################################################DELETE THE FOLDER OF THE STUDY##################################################
# for study in "${studies[@]}"; do
# rm -r "${path}""${study:?}/"* >/dev/null 
# mkdir -p "${path}""${study}"/nucl "${path}""${study}"/prot "${path}""${study}""${res_file}"/treefolder "${path}""${study}"/legend "${path}""${study}"/blast/database/ "${path}""${study}"/Align/
# done

# echo -e "Serogroup\tRegion\tOrder\nNa\t#979a9a" > "${path}"/All_legend_caption.txt
# # #Common scripts for all studies
# echo The sequences are currently being downloaded ...


# ######################################DOWNLOAD THE DATA IN GENBANK FORMAT FOR ALL THE ID FOR EACH STUDY##################################
# for i in "${!NUID[@]}"; do   
#     study="${studies[$i]}"
#     echo Download the nucleic sequences for the study : "${study}"
#     IFS=$';' read -ra IDprot <<< "${NUID[$i]}"
#     rm "${path}""${study}"/sequences_"${study}".gb

#    for ID in "${IDprot[@]}"; do 
#         ID=$(echo "$ID" | tr -d '[:space:]')
#         esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb
#         cat "${path}""${study}"/"$ID".gb >> "${path}""${study}"/sequences_"${study}".gb
#         rm "${path}""${study}"/"$ID".gb >/dev/null
#     done

# done
##########################################################SPECIFIQUE OUTGROUP###########################################################

# for i in "${!NUIDspeOG[@]}"; do   
#     study="${studiesspeOG[$i]}"
#     IFS=$';' read -ra IDprotspeOG <<< "${NUIDspeOG[$i]}"

#     for ID in "${IDprotspeOG[@]}"; do 

#         esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb

#     done

# done


############################################EXTRACT ALL THE DATA FROM THE GENBANK FILES#################################################
echo The captions, metadata files and files containing all the sequences are currently being created.

# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}/nucl/*" "${path}""${study:?}/prot/*" >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"genbank_extract.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${names[*]}" "${NUIDspeOG[*]}" "${email}" "${proteins[*]}"  "${proteinsspe_OG[*]}"   "${OG_prot[*]}" "${diff_protname[*]}" "${diff_protname_OG[*]}" "${order[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${Allstudies[*]}" "${allID[*]}" "${namespepro[*]}" "${namespeproOG[*]}"



# for i in "${!studies[@]}"; do

#     IFS=$';' read -ra prot <<< "${proteins[i]}"
#     IFS=$';' read -ra pos <<< "${position[i]}"
#     inc=0        
#     rm "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta "${path}""${studies[i]}"/"${target_list[i]}"_ref_DB
    
#     for p in "${!prot[@]}"; do
#         ss=${p}+${inc}
#         inc=${inc}+1
#         echo "${pos[ss]}"  "${pos[ss+1]}"
#         awk -v start="${pos[ss]}" -v end="${pos[ss+1]}" -v id="${target_list[i]}" -v proteine="${prot[p]}" '
#         !/^>/ { seq = seq $0 }          # concatène toutes les lignes de séquence
#         END {
#         print ">"  id "_" proteine 
#         print substr(seq, start, end - start + 1)
#         }
#         ' "${path}""${studies[i]}"/prot/"${target_list[i]}"_prot.fasta >> "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta
        
#     done
#         makeblastdb -in "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta \
#                     -title "${target_list[i]}"_ref.fasta \
#                     -dbtype prot \
#                     -out "${path}""${studies[i]}"/"${target_list[i]}"_ref_DB \
#                     -parse_seqids
# done

# for study in "${studies[@]}"; do
# mkdir -p "${path}""${studies[i]}"/Align
# rm -f "${path}""${study}"/Align/* >/dev/null 
# done

# echo "${path}""${path_script}"prot_position.py
# /usr/bin/python3 "${path}""${path_script}"prot_position.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${target_list[*]}" "${names[*]}" "${Allstudies[*]}"



# # # #####################################MESURE THE DIVERGENCE BETWEEN ALL THE SEQUENCES FOR EACH STUDY######################################
# ./"${path_script}"DNA_dist.sh "${path}" "${studies[*]}" 
# ######EQUIVALENT DNADIST - PROTDIST
# # # ##################################################DELETE THE SEQUENCES THAT ARE TOO CLOSE#################################################
# echo The duplicate sequences are being deleted
# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}"/nucl_sdbl_"${study}".fasta "${path}""${study:?}"/prot_sdbl_"${study}".fasta  >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "${studies[*]}" "${keep_prot[*]}" "${ths[*]}" "${target_list[*]}" 
# # # ######################################################SEPERATE THE PROTEINES##############################################################
# echo the proteins are being seperated
# for study in "${studies[@]}"; do
# echo 
# rm -f "${path}""${study:?}""${res_file}"/* >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"test_script_mark_spe.py "${path}" "${studies[*]}" "${position[*]}" "${names[*]}" "${positionspeOG[*]}" "${start[*]}" "${startspe_OG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${res_file}" "${studiesspeOG[*]}" "${NUIDspeOG[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}"



# # # # ##############################################ADD THE SPECIFIC OUTGROUP IN EACH FILE######################################################
# echo The outgroup of each studies are being added
# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}""${res_file}"/*_OG.fasta "${path}""${study:?}""${res_file}"/*_OG.aln-fasta  >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"merge_spe_prot.py "${path}" "${names[*]}" "${speprot[*]}" "${Type_seq[*]}" "${res_file}" "${studies[*]}" "${Allstudies[*]}" "${for_prot[*]}" 

# # # #####################################CREATE A TABLE WITH THE DIVERGENCE OF THE SPECIFIC PROTEINS############################################
# for study in "${studies[@]}"; do
#     mkdir -p "${path}""${study}"/Align_all
# done

# /usr/bin/python3 "${path}""${path_script}"protdist.py "${path}" "${studies[*]}" "${speprot[*]}" "${names[*]}" "${path_script}" "${res_file}" "${path_script}"


# # ##########################################"#LAUCH TRIMAL FOR EACH SPECIFIC PROTEINES OF EACH STUDY########################################
# speprotsse=$(echo "${speprot[@]}" | tr -d '" ' | cut -d'=' -f1 | paste -sd, -)
# IFS=$',' read -ra speprotsep <<< "${speprotsse[@]}"

# for i in "${!speprotsep[@]}"; do
#     if [[ -n "${speprotsep[$i]}" ]]; then
#         pairs+=("${speprotsep[$i]},${Allstudies[$i]}")
#     fi
# done

# # Supprimer les doublons
# mapfile -t unique_pairs < <(printf "%s\n" "${pairs[@]}" | sort -u)



# for pair in "${unique_pairs[@]}" ; do 
#     IFS=$',' read -ra protstud <<< "${pair[@]}"
#     if [[ " ${trim[*]} " =~ ${protstud[1]} ]]; then

#                 trimal  -in "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta" \
#                     -out "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" \
#                     -automated1 \
#                     -fasta 
#                 trimal  -in "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta" \
#                         -out "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.phy" \
#                         -automated1 \
#                         -phylip 
            

#             else 
#                 echo "convert" "${protstud[1]}"  "${protstud[0]}"
#                 seqkit grep -vs -r -p "^[Nn\-]+$" "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta">"${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta"
#                 seqret -sequence "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" -outseq "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.phy" -osformat2 phylip


#     fi

# done 




################################LAUCH TRIMAL FOR EACH PROTEINES OF EACH STUDY#############################################################
for study in "${!studies[@]}"; do

        IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
        for prot in "${sepprot[@]}"; do
            echo trimal pour la proteine "${prot}" pour l etude "${studies[$study]}"
            if [[ " ${trim[*]} " =~ ${studies[study]} ]]; then

                trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
                    -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" \
                    -automated1 \
                    -fasta 
                trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
                        -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" \
                        -automated1 \
                        -phylip 
            

            else 
                echo "convert" "${studies[$study]}"  "${prot}"
                seqkit grep -vs -r -p "^[Nn\-]+$" "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta">"${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta"
                seqret -sequence "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -outseq "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" -osformat2 phylip


             fi
        done
    
done


#################################################LAUCH IQTREE2 FOR EACH PROTEINES FOR EACH STUDY##########################################

for study in "${!studies[@]}"; do
    printf "The best model for each proteic sequence will be writen \n">"${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
    IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    for prot in "${sepprot[@]}"; do
            echo IQtree2 pour la proteine "${prot}" pour l étude "${studies[$study]}"
            { echo "${prot}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${studies[$study]}""${res_file}"treefolder/"${prot}"_"${studies[$study]}"_prot_trimal.phy ;}>> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
            iqtree2 -redo -s "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO -T AUTO |grep 'Best-fit model:' >> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
            printf "\n" >> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
    done
done


######################################LAUCH IQTREE2 FOR EACH  SPECIFIC PROTEINES FOR EACH STUDY###########################################
for pair in "${unique_pairs[@]}" ; do 
    IFS=$',' read -ra protstud <<< "${pair[@]}"

    echo IQtree2 pour la proteine "${protstud[0]}" pour l étude "${protstud[1]}"
            { echo "${protstud[0]}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${protstud[1]}""${res_file}"treefolder/"${protstud[0]}"_"${protstud[1]}"_prot_trimal.phy ;}>> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
            iqtree2 -redo -s "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO -T AUTO |grep 'Best-fit model:' >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
            printf "\n" >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt

done
Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "It is finally over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"

#snakemake
#nextflow