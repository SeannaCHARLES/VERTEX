#!/bin/bash
SECONDS=0

##########################################EXTRACT ALL THE INFORMATION FROM THE INFORMATION FILE##########################################
echo All the information are being extracted 

awk -F',' '$0 !~ /^#/ && $13 == "" { $13 = "none*nn=nn" } { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt

awk -F',' '$0 !~ /^#/ && $4 == "" { $4 = "N" } { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt
awk -F',' '$0 !~ /^#/ && $11 == "" { $11 = "0.05" } { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt
awk -F',' '$0 !~ /^#/ && $14 == "" { $14 = "Yes" } { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt

awk -F',' '$0 ~ /^#/ && NR==2 && $2=="" {$2="/"} { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt
awk -F',' '$0 ~ /^#/ && NR==3 && $2=="" {$2="/res_phylo/"} { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt
awk -F',' '$0 ~ /^#/ && NR==4 && $2=="" {$2="Phylogenetic/Script/"} { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt
awk -F',' '$0 ~ /^#/ && NR==5 && $2=="" {$2="nucl;prot"} { print }' OFS=',' all_information.txt > tmp && mv tmp all_information.txt


while awk -F',' 'NR < 6 { next }  NR == 6 {if ($0 ~ /^#/ && $2 == "") {exit 0} else {exit 1}}' OFS=',' "all_information.txt"; do

    echo "Please enter your email address :"
    read -p "email address : " INFORMATION

    if [ -z "$INFORMATION" ]; then
        echo "You cannot leave this field blank."

    else
        echo "saved email address: $INFORMATION"
        
        awk -F',' -v email="$INFORMATION" '$0 ~ /^#/ && NR==6 && $2=="" {$2=email} { print }' OFS=',' "all_information.txt" > "tmp_$$"
        mv "tmp_$$" "all_information.txt"
        

    fi
done

awk -F',' '{ $5 = tolower($5); $6 = tolower($6); $7 = tolower($7);$13 = toupper($13);$15 = tolower($15); OFS=","; print }' "all_information.txt" > tmp && mv tmp all_information.txt

while IFS=$',' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 _; do
    
    [[ $col1 =~ ^# ]] && continue

    names+=("$col2") #All the study and specific study name
    Allstudies+=("$col1") #All the Project name
    for_prot+=("$col5") # All the proteins for all the project
    
    
    if [[ $col4 == N* ]]; then
        studies+=("$col1")  #Name of the study that are not specific 
        
        mapfile -t L1 < <(awk  '{ for (i=1; i<=NF; i++) print $i }' "${col3}") 
        NUID+=("$(printf '%s ' "${L1[@]}")")

        target_list+=("$col8")  #References organism
        position+=("$col9") #References organism protein positions 
        start+=("$col10")   #Start of the protein sequence in the nucleic acid sequence
        ths+=("$col11") #Threshold for DNAdist 
        proteins+=("$col6") # Name of the proteins of all the organism
        namespepro+=("$col7")  # Name of the proteins of the
        if [[ $col14 == Y* ]]; then
            trim+=("$col1")
        fi
        blastapprox+=("$col15")
 
    else
        studiesspeOG+=("$col1") #Studies associé aux names 
           #Nom des 'sous unités' de chaque studies ex. spécifique outgroup
        #files_NUIDspeOG+=("$col3")    #Identifiants des outgroup spécifique
        mapfile -t L1 < <(awk  '{ for (i=1; i<=NF; i++) print $i }' "${col3}") 
        NUIDspeOG+=("$(printf '%s ' "${L1[@]}")")
        OG_prot+=("$col5")   #Proteines associées aux outgroup specifique
        positionspeOG+=("$col9")    #Positions des proteines de l'outgroup
        startspe_OG+=("$col10")  #Début de la séquence proteiques sur la séquence nucleique
        proteinsspe_OG+=("$col6")
     
        namespeproOG+=("$col7")
    fi
    if [[ -n $col7 ]]; then
        studiesspeprot+=("$col7")

        
    fi
done < all_information.txt



path=$(awk -F',' '$0 ~ /^#/ && NR==2 {print $2}' all_information.txt)

res_file=$(awk -F',' '$0 ~ /^#/ && NR==3 {print $2}' all_information.txt)

path_script=$(awk -F',' '$0 ~ /^#/ && NR==4 {print $2}' all_information.txt)

Type_seq=$(awk -F',' '$0 ~ /^#/ && NR==5 {print $2}' all_information.txt)

email=$(awk -F',' '$0 ~ /^#/ && NR==6 {print $2}' all_information.txt)

keep_prot=$(awk -F"," '$0 !~ /^#/ {print $12}' all_information.txt)

speprot=$(awk -F"," '$0 !~ /^#/ {print "\"",$7,"\""}' all_information.txt)

orgaWIWO=$(awk -F"," '$0 !~ /^#/ {gsub(/\xc2\xa0/, "", $13); print $13}' all_information.txt | awk -F"=" '{print $2}')

INOUT=$(awk -F"," '$0 !~ /^#/ {print $13}' all_information.txt|awk -F"=" '{print $1}'| awk -F "*" '{print $1}')

protINOUT=$(awk -F"," '$0 !~ /^#/{print $13}' all_information.txt|awk -F"=" '{print $1}'| awk -F "*" '{print $2}')


for element in "${speprot[@]}"; do
    if [ "$element" != " " ]; then
        speprot_empty=false
    
    fi
done

if $speprot_empty ; then
speprot=()
fi
allID=("${NUID[@]}" "${NUIDspeOG[@]}")

logfile () {

    rm -f "${path}""$2"/list_id_"$1".txt
    echo -e "After $1.py/sh : " >> "${path}""$2"/log_file_"$2".txt
    
    echo -e "Number of sequences : " >> "${path}""$2"/log_file_"$2".txt
    
    bioawk -c fastx '{print $name}' "$3" | wc -l >> "${path}""$2"/log_file_"$2".txt

    echo -e "length of the Alignement : " >> "${path}""$2"/log_file_"$2".txt
    bioawk -c fastx 'NR==1 {print length($seq)}' "$3" >> "${path}""$2"/log_file_"$2".txt

    echo -e "Deleted sequences : " >> "${path}""$2"/log_file_"$2".txt
    bioawk -c fastx '{print $name}' "$3" > "${path}""$2"/list_id_"$1".txt 
    sort "${path}""$2"/list_id_"$1".txt  > "${path}""$2"/list_id_"$1"_sort.txt 
    rm -f "${path}""$2"/list_id_"$1".txt


    diff  "${path}""$2"/list_id_"$4"_sort.txt "${path}""$2"/list_id_"$1"_sort.txt |grep "^<" >> "${path}""$2"/log_file_"$2".txt

    Hours=$(("${SECONDS}" / 3600))
    Minutes=$((("${SECONDS}" % 3600) / 60))
    Sec=$(("${SECONDS}" % 60))
    
    echo -e "-------------------------------------------------- ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds" >> "${path}""$2"/log_file_"$2".txt



}






################################################DELETE THE FOLDER OF THE STUDY#################################################

# for study in "${studies[@]}"; do

# rm -fr "${path}""${study:?}"/* >/dev/null 
# mkdir -p "${path}""${study}"/nucl "${path}""${study}"/prot "${path}""${study}""${res_file}"treefolder "${path}""${study}"/legend "${path}""${study}"/blast/database/ "${path}""${study}"/Align/ "${path}""${study}"/All_table/

# echo -e "legend\tcolor\nNa\t#979a9a" > "${path}""${study}"/all_legend_caption.txt
# done

# #Common scripts for all studies
# echo The sequences are currently being downloaded ...


###################################DOWNLOAD THE DATA IN GENBANK FORMAT FOR ALL THE ID FOR EACH STUDY####################################


# for i in "${!NUID[@]}"; do   
#     study="${studies[$i]}"
#     echo Download the nucleic sequences for the study : "${study}"
    
#     rm -f "${path}""${study}"/sequences_"${study}".gb >/dev/null 
#     read -ra IDprot <<< "${NUID[$i]}"

#    for ID in "${IDprot[@]}"; do 
#         esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb
#         sleep 2
#         cat "${path}""${study}"/"$ID".gb >> "${path}""${study}"/sequences_"${study}".gb
#         rm -f "${path}""${study}"/"$ID".gb >/dev/null
#     done

#     echo -e "Downloaded sequences  : " > "${path}""${study}"/log_file_"${study}".txt
#     awk '/^LOCUS/ {count++} END {print count}' "${path}""${study}"/sequences_"${study}".gb >> "${path}""${study}"/log_file_"${study}".txt
    
    
#     Hours=$(("${SECONDS}" / 3600))
#     Minutes=$((("${SECONDS}" % 3600) / 60))
#     Sec=$(("${SECONDS}" % 60))
#     echo -e "-------------------------------------------------- ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds" >> "${path}""${study}"/log_file_"${study}".txt

# done

######################################DOWNLOAD THE DATA IN GENBANK FORMAT FOR SPECIFIQUE OUTGROUP#######################################

# if [ "$speprot_empty" = "false" ]; then
# for i in "${!NUIDspeOG[@]}"; do   
#     study="${studiesspeOG[$i]}"
#     #IFS=$';' read -ra IDprotspeOG <<< "${NUIDspeOG[$i]}"
#     echo SPe ID download en cours
#     read -ra IDprotspeOG <<< "${NUIDspeOG[$i]}"
#     for ID in "${IDprotspeOG[@]}"; do 

#         esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb
#          sleep 2
#     done

# done
# fi

##############################################EXTRACT ALL THE DATA FROM THE GENBANK FILES###############################################

# echo The captions, metadata files and files containing all the sequences are currently being created.

# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}/nucl/*" "${path}""${study:?}/prot/*" >/dev/null 
# done


# /usr/bin/python3 "${path}""${path_script}"genbank_extract.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${names[*]}" "${NUIDspeOG[*]}" "${email}" "${proteins[*]}"  "${proteinsspe_OG[*]}"   "${OG_prot[*]}"  "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${Allstudies[*]}" "${allID[*]}" "${namespepro[*]}" "${namespeproOG[*]}"

# for study in "${studies[@]}"; do
#     echo -e "After genbank_extract.py : " >> "${path}""${study}"/log_file_"${study}".txt
#     echo -e "Number of sequences : " >> "${path}""${study}"/log_file_"${study}".txt

#     bioawk -c fastx '{print $name}' "${path}""${study}"/prot_all_align_"${study}".aln-fasta | wc -l >> "${path}""${study}"/log_file_"${study}".txt

#     bioawk -c fastx '{print $name}' "${path}""${study}"/prot_all_align_"${study}".aln-fasta >  "${path}""${study}"/list_id_genbank.txt 

#     sort "${path}""${study}"/list_id_genbank.txt  > "${path}""${study}"/list_id_genbank_sort.txt
#     rm -f "${path}""${study}"/list_id_genbank.txt

#     echo -e "length of the Alignement : " >> "${path}""${study}"/log_file_"${study}".txt
#     bioawk -c fastx 'NR==1 {print length($seq)}' "${path}""${study}"/prot_all_align_"${study}".aln-fasta >> "${path}""${study}"/log_file_"${study}".txt

#     Hours=$(("${SECONDS}" / 3600))
#     Minutes=$((("${SECONDS}" % 3600) / 60))
#     Sec=$(("${SECONDS}" % 60))
#     echo -e "-------------------------------------------------- ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds" >> "${path}""${study}"/log_file_"${study}".txt

# done



########################################COMPLETE THE MISSING PROTEIN POSTIONS##########################################

# for i in "${!studies[@]}"; do
#     IFS=$';' read -ra prot <<< "${proteins[i]}"
#     IFS=$';' read -ra pos <<< "${position[i]}"
#     inc=0        
#     rm -f "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta "${path}""${studies[i]}"/"${target_list[i]}"_ref_DB >/dev/null
#     echo "${path}""${studies[i]}"/prot/"${target_list[i]}"_prot.fasta
#     echo "${proteins[i]}"
#     echo "${position[i]}"


#     for p in "${!prot[@]}"; do
#         ss=${p}+${inc}
#         inc=${inc}+1
#         echo "${pos[ss]}"  "${pos[ss+1]}"
#         echo "${target_list[i]}" 
#         echo "${prot[p]}"
#         awk -v start="${pos[ss]}" -v end="${pos[ss+1]}" -v id="${target_list[i]}" -v proteine="${prot[p]}" '
#         !/^>/ { seq = seq $0 }         
#         END {
#         print ">"  id "_" proteine 
#         print substr(seq, start, end - start + 1)
#         }' "${path}""${studies[i]}"/prot/"${target_list[i]}"_prot.fasta >> "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta
        
#     done
#         makeblastdb -in "${path}""${studies[i]}"/"${target_list[i]}"_ref.fasta \
#                     -title "${target_list[i]}"_ref.fasta \
#                     -dbtype prot \
#                     -out "${path}""${studies[i]}"/"${target_list[i]}"_ref_DB \
#                     -parse_seqids
#     echo 
# done

# for study in "${studies[@]}"; do

#         mkdir -p "${path}""${study}"/Align
#         rm -f "${path}""${study}"/Align/* >/dev/null 
#         echo "${path}""${path_script}"prot_position.py
#         /usr/bin/python3 "${path}""${path_script}"prot_position.py "${path}" "${studies[*]}" "${studiesspeOG[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}" "${target_list[*]}" "${names[*]}" "${Allstudies[*]}" "${blastapprox[*]}"
    

# done 





#####################################MESURE THE DIVERGENCE BETWEEN ALL THE SEQUENCES FOR EACH STUDY######################################

# ./"${path_script}"DNA_dist.sh "${path}" "${studies[*]}" 

################################################DELETE THE SEQUENCES THAT ARE TOO CLOSE#################################################

# echo The duplicate sequences are being deleted
# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}"/nucl_sdbl_"${study}".fasta "${path}""${study:?}"/prot_sdbl_"${study}".fasta  >/dev/null 
# echo -e "After delete_dbl.py : " >> "${path}""${study}"/log_file_"${study}".txt
# echo -e "Deleted sequences : " >> "${path}""${study}"/log_file_"${study}".txt
# done

/usr/bin/python3 "${path}""${path_script}"delete_dbl.py "${path}" "${studies[*]}" "${keep_prot[*]}" "${ths[*]}" "${target_list[*]}" 
for study in "${studies[@]}"; do

    logfile "delete_dbl" "${study}" "${path}""${study}"/prot_sdbl_"${study}".fasta "genbank"
done

################################################SEPERATE THE PROTEINES#################################################


# echo the proteins are being seperated

# for study in "${studies[@]}"; do
# # IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
# # echo "${sepprot[*]}"
# rm -f "${path}""${study:?}""${res_file}"* >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"sep_prot.py "${path}" "${studies[*]}" "${position[*]}" "${names[*]}" "${positionspeOG[*]}" "${start[*]}" "${startspe_OG[*]}" "${proteins[*]}" "${proteinsspe_OG[*]}" "${res_file}" "${studiesspeOG[*]}" "${NUIDspeOG[*]}" "${orgaWIWO[*]}" "${INOUT[*]}" "${protINOUT[*]}"


# for study in "${studies[@]}"; do
#     #echo -e "After sep_prot.py : " >> "${path}""${study}"/log_file_"${study}".txt
    
#     for prot in "${sepprot[@]}"; do
#         echo -e "${prot} : " >> "${path}""${study}"/log_file_"${study}".txt
#         logfile  "sep_prot" "${study}" "${path}""${study}""${res_file}""${prot}"_"${study}"_prot.fasta "delete_dbl"
#     done
# done


#################################################DON'T COMMENT THAT PART#################################################

for study in "${!studies[@]}"; do
    
    IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    maxprot=${#sepprot[*]}
    for (( i=3; i<=(maxprot*2)+2; i+=2 )); do
        resultat=$(awk -F"\t" -v col="$i" 'NR>1 && ($col!="nan" || $(col+1)!="nan"){found=1;exit} END{if(found)print "not only nan";else print "only nan"}' "${path}${studies[$study]}/pos_prot_final_${studies[$study]}.txt")
        todelprot=$(head -n1 "${path}${studies[$study]}/pos_prot_final_${studies[$study]}.txt" | awk -v col="$i" -F"\t" '{print $col}' | sed 's/Start//g')

        if [ "$resultat" = "only nan" ]; then
            
            echo "The protein ${todelprot} does not have any positions and therefor will be deleted for the study ${studies[$study]}."
            echo "The protein ${todelprot} does not have any positions and therefor will be deleted for the study ${studies[$study]}." >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
            newprotlist=""
            #newprotlistsep=""
            for element in "${sepprot[@]}"; do
                if [[ "$element" != "$todelprot" ]]; then
                    if [[ -n "$newprotlist" ]]; then
                        newprotlist+=";"
                        #newprotlistsep+=" "
                    fi
                    newprotlist+="$element"
                    #newprotlistsep+="$element"
                fi
            done
            proteins["$study"]="${newprotlist[*]}"
            #sepprot[*]="${newprotlistsep[*]}"

            newprotlist=""
            IFS=$';' read -ra sepprotspe <<< "${proteinsspe_OG[$study]}"
            for element in "${sepprotspe[@]}"; do
                if [[ "$element" != "$todelprot" ]]; then
                    if [[ -n "$newprotlist" ]]; then
                        newprotlist+=";"
                    fi

                    newprotlist+="$element"
                fi
            done
            proteinsspe_OG["$study"]="${newprotlist[*]}"


            newprotlist=""
            
            IFS=$';' read -ra sepprotfor <<< "${for_prot[$study]}"
            #newprotlistsep=""
            for element in "${sepprotfor[@]}"; do
                if [[ "$element" != "$todelprot" ]]; then
                    if [[ -n "$newprotlist" ]]; then
                        newprotlist+=";"
                        #newprotlistsep+=" "
                    fi
                    newprotlist+="$element"
                    #newprotlistsep+="$element"
                fi
            done
            for_prot["$study"]="${newprotlist[*]}"
            
        fi
    done

done



for study in "${!studies[@]}"; do

IFS=$';' read -ra sepprot <<< "${proteins[$study]}" 

    for prot in "${sepprot[@]}"; do
        echo -e "${prot} : " >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
        logfile  "sep_prot" "${studies[$study]}" "${path}""${studies[$study]}""${res_file}""${prot}"_"${studies[$study]}"_prot.fasta "delete_dbl"
    done
done

#################################################ADD THE SPECIFIC OUTGROUP IN EACH FILE#################################################

# echo The outgroup of each studies are being added
# for study in "${studies[@]}"; do
# rm -f "${path}""${study:?}""${res_file}"*_OG.fasta "${path}""${study:?}""${res_file}"*_OG.aln-fasta  >/dev/null 
# done

# /usr/bin/python3 "${path}""${path_script}"merge_spe_prot.py "${path}" "${names[*]}" "${speprot[*]}" "${Type_seq[*]}" "${res_file}" "${studies[*]}" "${Allstudies[*]}" "${for_prot[*]}" 

for study in "${!studies[@]}"; do
    
   IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    
    for prot in "${sepprot[@]}"; do
        echo -e "${prot} : " >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
        logfile  "merge_spe_prot" "${studies[$study]}" "${path}""${studies[$study]}""${res_file}""${prot}"_"${studies[$study]}"_prot_OG.aln-fasta "sep_prot"
    done
done
    



##########################################LAUCH TRIMAL FOR EACH SPECIFIC PROTEINES OF EACH STUDY##########################################

# if [ "${speprot[0]}" != '"  "' ];then
#     speprotsse=$(echo "${speprot[@]}" | tr -d '" ' | cut -d'=' -f1 | paste -sd, -)
#     IFS=$',' read -ra speprotsep <<< "${speprotsse[@]}"

#     for i in "${!speprotsep[@]}"; do
#         if [[ -n "${speprotsep[$i]}" ]]; then
#             pairs+=("${speprotsep[$i]},${Allstudies[$i]}")
#         fi
#     done


#     mapfile -t unique_pairs < <(printf "%s\n" "${pairs[@]}" | sort -u)



#     for pair in "${unique_pairs[@]}" ; do 
#         IFS=$',' read -ra protstud <<< "${pair[@]}"
#         if [[ " ${trim[*]} " =~ ${protstud[1]} ]]; then

#                     trimal  -in "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta" \
#                         -out "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" \
#                         -automated1 \
#                         -fasta 
#                     trimal  -in "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta" \
#                             -out "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.phy" \
#                             -automated1 \
#                             -phylip 
            

#                 else 
#                     echo "convert" "${protstud[1]}"  "${protstud[0]}"
#                     seqkit grep -vs -r -p "^[Nn\-]+$" "${path}${protstud[1]}${res_file}${protstud[0]}_${protstud[1]}_prot_OG.aln-fasta">"${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta"
#                     seqret -sequence "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" -outseq "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.phy" -osformat2 phylip


#         fi

#     done 
# fi



##############################################LAUCH TRIMAL FOR EACH PROTEINES OF EACH STUDY##############################################

# for study in "${!studies[@]}"; do
#   IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
        
#         for prot in "${sepprot[@]}"; do
#             echo trimal pour la proteine "${prot}" pour l etude "${studies[$study]}"
#             if [[ " ${trim[*]} " =~ ${studies[study]} ]]; then

#                 trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
#                     -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" \
#                     -automated1 \
#                     -fasta \
#                     -htmlout "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal_report.html"

#                 trimal  -in "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta" \
#                         -out "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" \
#                         -automated1 \
#                         -phylip \
#                         -htmlout "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal_report_phy.html"
            

#             else 
#                 echo "convert" "${studies[$study]}"  "${prot}"

#                 seqkit grep -vs -r -p "^[Nn\-]+$" "${path}${studies[$study]}${res_file}${prot}_${studies[$study]}_prot_OG.aln-fasta">"${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta"

#                 seqret -sequence "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -outseq "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.phy" -osformat2 phylip


#              fi
#         done
    
# done

for study in "${!studies[@]}"; do
  IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    
    for prot in "${sepprot[@]}"; do
        logfile  "trimal" "${studies[$study]}" "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" "merge_spe_prot"
    done
done



#############################################LAUCH IQTREE2 FOR EACH PROTEINES FOR EACH STUDY#############################################

for study in "${!studies[@]}"; do
    printf "The best model for each proteic Alignements will be writen \n">>"${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
    IFS=$';' read -ra sepprot <<< "${proteins[$study]}"
    for prot in "${sepprot[@]}"; do
            echo IQtree2 pour la proteine "${prot}" pour l étude "${studies[$study]}"
            
            { echo "${prot}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Alignements size:", $2}' "${path}""${studies[$study]}""${res_file}"treefolder/"${prot}"_"${studies[$study]}"_prot_trimal.phy ; } >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt

            iqtree2 --redo -s "${path}${studies[$study]}${res_file}treefolder/${prot}_${studies[$study]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO -T AUTO |grep 'Best-fit model:' >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
            printf "\n" >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
echo -e "-------------------------------------------------- ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds" >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
    done
done


########################################LAUCH IQTREE2 FOR EACH  SPECIFIC PROTEINES FOR EACH STUDY########################################

if [ "$speprot_empty" = "false" ]; then
    for pair in "${unique_pairs[@]}" ; do 
        IFS=$',' read -ra protstud <<< "${pair[@]}"

        echo IQtree2 pour la proteine "${protstud[0]}" pour l étude "${protstud[1]}"
                { echo "${protstud[0]}" : ; awk 'NR==1 {print "Number of sequences:", $1, "Alignements size:", $2}' "${path}""${protstud[1]}""${res_file}"treefolder/"${protstud[0]}"_"${protstud[1]}"_prot_trimal.phy ; } >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
                iqtree2 --redo -s "${path}${protstud[1]}${res_file}treefolder/${protstud[0]}_${protstud[1]}_prot_trimal.fasta" -m MFP -bb 1000 -alrt 1000 -nt AUTO -T AUTO |grep 'Best-fit model:' >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
                printf "\n" >> "${path}""${protstud[1]}""${res_file}"Best_model_"${protstud[1]}".txt
echo -e "-------------------------------------------------- ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds" >> "${path}""${studies[$study]}"/log_file_"${studies[$study]}".txt
    done
fi
Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "End of IQtree2 you can already use the phylogenetic tree in newick format and the associated caption. "
#######################################CREATE A TABLE WITH THE DIVERGENCE OF THE SPECIFIC PROTEINS#######################################
echo "Computes a protein distance matrix for each protein and specific protein using Protdist this step migh take a while "
if [ "$speprot_empty" = "false" ]; then
for study in "${studies[@]}"; do
    mkdir -p "${path}""${study}"/Align_all
done

/usr/bin/python3 "${path}""${path_script}"protdist.py "${path}" "${studies[*]}" "${speprot[*]}" "${names[*]}" "${path_script}" "${res_file}" "${path_script}" "${proteins[*]}"  "${studiesspeprot[*]}"
fi
echo "It is over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"