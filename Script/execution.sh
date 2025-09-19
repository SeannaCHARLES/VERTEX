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
        namespepro+=("$col7")
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
rm -r "${path}""${study:?}/"* >/dev/null 
mkdir -p "${path}""${study}"/nucl "${path}""${study}"/prot "${path}""${study}""${res_file}"/treefolder "${path}""${study}"/legend "${path}""${study}"/blast/database/
done
echo -e "Serogroup\tRegion\tOrder\nNa\t#979a9a" > "${path}"/All_legend_caption.txt
# #Common scripts for all studies
 echo The sequences are currently being downloaded ...



for i in "${!NUID[@]}"; do   
    study="${studies[$i]}"
    echo Download the nucleic sequences for the study : "${study}"
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
for study in "${studies[@]}"; do
rm -f "${path}""${study:?}/nucl/*" "${path}""${study:?}/prot/*" >/dev/null 
done
/usr/bin/python3 "${path}""${path_script}"genbank_extract.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${names[*]}" "${NUIDspeOG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${email}"

./"${path_script}"DNA_dist.sh "${path}" "${studies[*]}" 

echo The duplicate sequences are being deleted
for study in "${studies[@]}"; do
rm -f "${path}""${study:?}"/nucl_sdbl_"${study}".fasta "${path}""${study:?}"/prot_sdbl_"${study}".fasta  >/dev/null 
done
/usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "${studies[*]}" "${keep_prot[*]}" "${ths[*]}" "${target_list[*]}" 


echo the proteins are being seperated
for study in "${studies[@]}"; do
rm -f "${path}""${study:?}""${res_file}"* >/dev/null 
done
/usr/bin/python3 "${path}""${path_script}"test_script_mark_spe.py "${path}" "${studies[*]}" "${position[*]}" "${names[*]}" "${positionspeOG[*]}" "${start[*]}" "${startspe_OG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${res_file}" "${studiesspeOG[*]}" "${NUIDspeOG[*]}"

#Blast pour avoir le % de similarité table pour chaque study
./"${path_script}"Blast_prot.sh "${path}" "${res_file}" "${studies[*]}" 

echo Lauch of blast
for study in "${!studies[@]}"; do
    IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    for prot in "${sepprot[@]}"; do        
            makeblastdb -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot.fasta"\
                        -title "${prot}_${studies[$study]}${prot}"\
                        -dbtype prot \
                        -out "${path}${studies[$study]}/blast/database/${prot}_DB"\
                        -parse_seqids
            bioawk -c fastx '{ print $name, length($seq) }' "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot.fasta" > "${path}${studies[$study]}/blast/${prot}_${studies[$study]}_prot_size.fasta"
    done
        prot=$(echo "${namespepro[$study]}" | cut -d'=' -f1 | tr -d ' ')
        echo "${prot}"
    for prot1 in "${sepprot[@]}"; do
        for prot2 in "${sepprot[@]}"; do
            blastp  -query "${path}${studies[$study]}${res_file}${prot1}_${studies[$study]}_prot.fasta" \
                    -db "${path}${studies[$study]}/blast/database/${prot2}_DB" \
                    -out "${path}${studies[$study]}/blast/${prot1}_${prot2}_blast.txt" \
                    -outfmt "6"
            awk '                                                                                         ARGIND==1 { qsize[$1]=$2; next }        # 1er fichier : tailles des query
                ARGIND==2 { ssize[$1]=$2; next }        # 2e fichier : tailles des subject
                ARGIND==3 {
                qs = ($1 in qsize) ? qsize[$1] : "0"
                ss = ($2 in ssize) ? ssize[$2] : "0"
                print $0, qs, ss
                                                }
            ' "${path}${studies[$study]}/blast/${prot1}_${studies[$study]}_prot_size.fasta"   "${path}${studies[$study]}/blast/${prot2}_${studies[$study]}_prot_size.fasta"   "${path}${studies[$study]}/blast/${prot1}_${prot2}_blast.txt" > "${path}${studies[$study]}/blast/${prot1}_${prot2}_blast_sizes.txt"

            awk '{print $0 "\t" (($10-$9+1)/$13*100) "\t" (($8-$7+1)/$14*100)}'   "${path}${studies[$study]}/blast/${prot1}_${prot2}_blast_sizes.txt" | sed '1iQseqid\tSseqid\tpident\tlength\tmismatch\tgapopen\tQstart\tQend\tSstart\tSend\tevalue\tbitscore\tQlength\tSlength\tQcoverage\tScoverage'> "${path}${studies[$study]}/blast/${prot1}_${prot2}_blast_table.txt"

        done
    done

done

echo The proteines that are presented or not in certains proteins are being treated
/usr/bin/python3 "${path}""${path_script}"WIWO_protein.py "${path}" "${names[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${res_file}" "${Type_seq[*]}"

# #merge = OG
# echo The outgroup of each studies are being added
# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}"/"${res_file}"/*_OG.fasta "${path}""${study:?}"/"${res_file}"/*_OG.aln-fasta  >/dev/null 
# done
# /usr/bin/python3 "${path}""${path_script}"merge_spe_prot.py "${path}" "${names[*]}" "${speprot[*]}" "${Type_seq[*]}" "${res_file}" "${studies[*]}" "${Allstudies[*]}" "${for_prot[*]}" 
# speprotsse=$(echo "${speprot[@]}" | tr -d '" ' | cut -d'=' -f1 | paste -sd, -)

# IFS=$',' read -ra speprotsep <<< "${speprotsse[@]}"
# for i in "${!speprotsep[@]}"; do
#     if [[ -n "${speprotsep[$i]}" ]]; then
#         pairs+=("${speprotsep[$i]},${Allstudies[$i]}")
#     fi
# done

# # Supprimer les doublons
# mapfile -t unique_pairs < <(printf "%s\n" "${pairs[@]}" | sort -u)


# for study in "${!studies[@]}"; do

#         IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
#         for prot in "${sepprot[@]}"; do
#             echo trimal pour la proteine "${prot}" pour l etude "${studies[$study]}"
#             if [[ " ${trim[*]} " =~ ${studies[study]} ]]; then

#                 trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
#                     -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" \
#                     -automated1 \
#                     -fasta 
#                 trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
#                         -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" \
#                         -automated1 \
#                         -phylip 
            

#             else 
#                 echo "convert" "${studies[$study]}"  "${prot}"
#                 seqkit grep -vs -r -p "^[Nn\-]+$" "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta">"${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta"
#                 seqret -sequence "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -outseq "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" -osformat2 phylip


#              fi
#         done
    
# done

# for pair in "${unique_pairs[@]}" ; do 
# IFS=$',' read -ra protstud <<< "${pair[@]}"
# if [[ " ${trim[*]} " =~ ${protstud[1]} ]]; then

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


#              fi

# done 

# for study in "${!studies[@]}"; do
# printf "The best model for each proteic sequence will be write \n">"${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
#     IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
#     for prot in "${sepprot[@]}"; do
#             echo IQtree2 pour la proteine "${prot}" pour l étude "${studies[$study]}"
#             { echo "${prot}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${studies[$study]}""${res_file}"treefolder/"${prot}"_"${studies[$study]}"_prot_trimal.phy ;}>> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
#             iqtree2 -redo -s "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
#             printf "\n" >> "${path}""${studies[$study]}""${res_file}"Best_model_"${studies[$study]}".txt
#     done
# done
# for pair in "${unique_pairs[@]}" ; do 
# IFS=$',' read -ra protstud <<< "${pair[@]}"
# echo IQtree2 pour la proteine "${protstud[0]}" pour l étude "${protstud[1]}"
#             { echo "${protstud[0]}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Sequences size:", $2}' "${path}""${protstud[1]}""${res_file}"treefolder/"${protstud[0]}"_"${protstud[1]}"_prot_trimal.phy ;}>> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
#             iqtree2 -redo -s "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
#             printf "\n" >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt

# done

# Hours=$(("${SECONDS}" / 3600))
# Minutes=$((("${SECONDS}" % 3600) / 60))
# Sec=$(("${SECONDS}" % 60))

# echo "It is finally over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"