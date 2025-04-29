#!/bin/bash

SECONDS=0
path="/home/seannacharles/Documents/"
studies=("Peribunyaviridae" "CCHFV" "Nairoviridae")
Type_seq=("prot" "nucl")
proteins=("NSm" "Gn" "Gc" "NSmGn")

for study in "${studies[@]}"; do
mkdir -p "${path}""${study}"/M_nucl "${path}""${study}"/M_prot "${path}""${study}"/seg_M_phylo
done
#Common scripts for all studies
echo The sequences are currently being downloaded ...
./Phylogenetic/Script/download_sequences.sh

rm "${path}"Metadata_studies.txt >/dev/null
echo -e "Virus name\tID itol tree\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteique sequence\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tStart Gn\tEnd Gn\tStart NSm\tEnd NSm\tStart Gc\tEnd Gc\tNucleique accession number\tSize Nucelique sequence\tStrand\tStrain\tGeographical region\tHost\tHost source\tReferences\n" >"${path}"Metadata_studies.txt

echo The captions, metadata files and files containing all the sequences are currently being created.
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/genbank_extract.py


for study in "${studies[@]}"; do
    echo Start treatment for"${study}" : 

    echo The position of the proteins of the "${study}" are being marked ...
    /usr/bin/python3 /home/seannacharles/Phylogenetic/Script/mark_prot_"${study}".py

    echo The proteins are being splited for the "${study}" ...
    /usr/bin/python3 /home/seannacharles/Phylogenetic/Script/split_protein_"${study}".py

    echo The NSm and Gn protein identifiers are currently being modified to avoid duplicate identifiers for the "${study}"...
    if [ "$study" == "Peribunyaviridae" ] ;then 
        { sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprotSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotSH.fasta ; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
        { sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnuclSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclSH.fasta ;} > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta

    else
        { sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnprot.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
        { sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnnucl.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta

    fi
    
    mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot_align.aln-fasta
    mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl_align.aln-fasta

    echo Non-informative regions are being removed and sequences that are too close will be deleted for the "${study}"...
    ./Phylogenetic/Script/trimal_"${study}".sh

    echo Identification of the best phylogenetic model is underway for the "${study}"
    echo The phylogenetic trees are also being created. This operation may take a few hours depending on the number of sequences. 

   for type in "${Type_seq[@]}" ; do
        echo for the "${type}" sequences...
        for name in "${proteins[@]}";do
            echo "${name}""${type}" : >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}"/seg_M_phylo/"${name}""${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt
        done
        if [ "$study" != "Peribunyaviridae" ] ;then 
            echo MLD"${type}": >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}"/seg_M_phylo/MLD"${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt

            echo GP38"${type}": >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt
            iqtree2 -s "${path}""${study}"/seg_M_phylo/GP38"${type}"_sdbl.aln-fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_"${study}".txt
        fi

done
    echo End of treatment for"${study}" :
done


Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "It is finally over excution time : ${Hours} Hours ${Minutes} Minutes ${Sec} Seconds"