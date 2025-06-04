#!/bin/bash
SECONDS=0

path="/home/seannacharles/"
path_script="Phylogenetic/Script/"
res_file="/res_phylo/"
studies=("Peribunyaviridae" "CCHFV" "Nairoviridae")
Type_seq=("prot" "nucl")
proteins=("NSm" "Gn" "Gc" "NSmGn")

for study in "${studies[@]}"; do
rm -r "${path}""${study:?}/"* 
mkdir -p "${path}""${study}"/nucl "${path}""${study}"/prot "${path}""${study}""${res_file}"
done

#Common scripts for all studies
 echo The sequences are currently being downloaded ...
 ./"${path_script}"download_sequences.sh "${path}" 


echo The captions, metadata files and files containing all the sequences are currently being created.
/usr/bin/python3 "${path}""${path_script}"genbank_extract.py "${path}"


for study in "${studies[@]}"; do
    echo Start treatment for "${study}" : 

    echo The position of the proteins of the "${study}" are being marked ...
    /usr/bin/python3 "${path}""${path_script}"mark_prot_"${study}".py "${path}" "${study}"

    echo The proteins are being splited for the "${study}" ...
    /usr/bin/python3 "${path}""${path_script}"split_protein_"${study}".py  "${path}" "${study}" "${res_file}"

    echo The NSm and Gn protein identifiers are currently being modified to avoid duplicate identifiers for the "${study}"...
    if [ "$study" == "Peribunyaviridae" ] ;then 
        { sed 's/^>/>N/' "${path}""${study}""${res_file}"NSmprotSH.fasta; sed 's/^>/>G/' "${path}""${study}""${res_file}"GnprotKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}""${res_file}"GnprotSH.fasta ; } > "${path}""${study}""${res_file}"NSmGnprot.fasta
        { sed 's/^>/>N/' "${path}""${study}""${res_file}"NSmnuclSH.fasta; sed 's/^>/>G/' "${path}""${study}""${res_file}"GnnuclKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}""${res_file}"GnnuclSH.fasta ;} > "${path}""${study}""${res_file}"NSmGnnucl.fasta

    else
        { sed 's/^>/>N/' "${path}""${study}""${res_file}"NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}""${res_file}"Gnprot.fasta; } > "${path}""${study}""${res_file}"NSmGnprot.fasta
        { sed 's/^>/>N/' "${path}""${study}""${res_file}"NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}""${res_file}"Gnnucl.fasta; } > "${path}""${study}""${res_file}"NSmGnnucl.fasta

    fi
    
    mafft --auto "${path}""${study}""${res_file}"NSmGnprot.fasta > "${path}""${study}""${res_file}"NSmGnprot_align.aln-fasta
    mafft --auto "${path}""${study}""${res_file}"NSmGnnucl.fasta > "${path}""${study}""${res_file}"NSmGnnucl_align.aln-fasta

    echo Non-informative regions are being removed and sequences that are too close will be deleted for the "${study}"...
    ./"${path_script}"trimal_"${study}".sh "${study}" "${path}" "${res_file}" "${path_script}"

    echo Identification of the best phylogenetic model is underway for the "${study}"
    echo The phylogenetic trees are also being created. This operation may take a few hours depending on the number of sequences. 

   for type in "${Type_seq[@]}" ; do
        echo for the "${type}" sequences...
        for name in "${proteins[@]}";do
            mafft --auto "${path}""${study}""${res_file}""${name}""${type}"_sdbl.fasta > "${path}""${study}""${res_file}""${name}""${type}"_sdbl.aln-fasta
            seqret -sequence "${path}""${study}""${res_file}""${name}""${type}"_sdbl.aln-fasta -outseq "${path}""${study}""${res_file}""${name}""${type}"_Final.phy -osformat2 phylip
            { echo "${name}""${type}" : ; echo awk 'NR==1' "${path}""${study}""${res_file}""${name}""${type}"_Final.phy ;}>> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}""${res_file}""${name}""${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
        done
        if [ "$study" != "Peribunyaviridae" ] ;then 
            mafft --auto "${path}""${study}""${res_file}"MLD"${type}"_sdbl.fasta > "${path}""${study}""${res_file}"MLD"${type}"_sdbl.aln-fasta
            seqret -sequence "${path}""${study}""${res_file}"MLD"${type}"_sdbl.aln-fasta -outseq "${path}""${study}""${res_file}"MLD"${type}"_Final.phy -osformat2 phylip
            { echo MLD"${type}": ; echo awk 'NR==1' "${path}""${study}""${res_file}""${name}""${type}"_Final.phy ; } >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}""${res_file}"MLD"${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            
            mafft --auto "${path}""${study}""${res_file}"GP38"${type}"_sdbl.fasta > "${path}""${study}""${res_file}"GP38"${type}"_sdbl.aln-fasta
            seqret -sequence "${path}""${study}""${res_file}"GP38"${type}"_sdbl.aln-fasta -outseq "${path}""${study}""${res_file}"GP38"${type}"_Final.phy -osformat2 phylip
            { echo GP38"${type}": ; echo awk 'NR==1' "${path}""${study}""${res_file}"GP38"${type}"_Final.phy ; } >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}""${res_file}"GP38"${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}""${res_file}"Best_model_"${study}".txt
        fi

done
    echo End of treatment for "${study}" :
done

#./Phylogenetic/Script/execution.sh: ligne 77: erreur de syntaxe près du symbole inattendu « ) »
#./Phylogenetic/Script/execution.sh: ligne 77: `/ 3600))'
Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "It is finally over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"