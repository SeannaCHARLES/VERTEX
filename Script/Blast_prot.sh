#!/bin/bash

path="$1"
res_file="$2"
_studies="$3"
IFS=' ' read -ra studies <<< "$_studies"
echo "${studies[*]}"
for study in "${studies[@]}"; do
    echo "${study} blast"

    # Étape 1 : Créer toutes les bases de données
    for fichier in "${path}${study}${res_file}"*prot_OG.fasta; do
        nom_fichier=$(basename "$fichier")
        start_file_name="${nom_fichier%prot_OG.fasta}"

        makeblastdb -in "${path}${study}${res_file}${start_file_name}prot_OG.fasta" \
                    -title "${start_file_name}" \
                    -dbtype prot \
                    -out "${path}${study}/blast/database/${start_file_name}DB" \
                    -parse_seqids
    done

    # Étape 2 : Obtenir la liste des fichiers de base de données
    mapfile -t db_files < <(find "${path}${study}/blast/database" -name "*_DB")

    # Étape 3 : Comparer toutes les paires uniques
    for ((i = 0; i < ${#db_files[@]}; i++)); do
        query_file="${db_files[i]}"
        query_name=$(basename "${query_file}")
        query_prefix="${query_name%_DB}"

        for ((j = i + 1; j < ${#db_files[@]}; j++)); do
            db_file="${db_files[j]}"
            db_name=$(basename "${db_file}")
            db_prefix="${db_name%_DB}"

            # Lancer les deux alignements : A vs B et B vs A si besoin, ou seulement un selon ton besoin
            # Ici on fait seulement A vs B pour éviter les doublons
            blastp -query "${path}${study}${res_file}${query_prefix}prot_OG.fasta" \
                   -db "${db_file}" \
                   -out "${path}${study}/blast/${query_prefix}_vs_${db_prefix}_blast.txt" \
                   -outfmt 7
        done
    done
done
