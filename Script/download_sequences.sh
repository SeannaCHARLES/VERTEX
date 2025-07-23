#!/bin/bash

path="$1"
allstudies="$2"
NUID="$3"
studiesspeOG="$4"
NUIDspeOG="$5"
IFS=' ' read -ra studies <<< "$allstudies"


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
    echo "${study}" specific OG
    IFS=$';' read -ra IDprotspeOG <<< "${NUIDspeOG[$i]}"
    for ID in "${IDprotspeOG[@]}"; do 
         esearch -db nucleotide -query "$ID" |efetch -format gb > "${path}""${study}"/"$ID".gb

    done
done

