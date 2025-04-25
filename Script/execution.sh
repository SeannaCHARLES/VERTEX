#!/bin/bash
SECONDS=0
path="/home/seannacharles/Documents/"

#Common scripts for all studies
echo The sequences are currently being downloaded ...
./Phylogenetic/Script/download_sequences.sh

echo The captions, metadata files and files containing all the sequences are currently being created.
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/genbank_extract.py

study="Perybunyaviridae"
echo Start of treatment for"${study}" : 

echo The position of the proteins of the "${study}" are being marked ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/mark_prot.py

echo The proteins are being splited for the "${study}" ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/split_protein_peri.py

echo The NSm and Gn protein identifiers are currently being modified to avoid duplicate identifiers for the "${study}"...
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprotSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotSH.fasta ; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnuclSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclSH.fasta ;} > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta

mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot_align.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl_align.aln-fasta

echo Non-informative regions are being removed and sequences that are too close will be deleted for the "${study}"...
./Phylogenetic/Script/trimal_peri.sh

echo Identification of the best phylogenetic model is underway for the "${study}"

echo For the proteique sequences...
echo this operation may take up to ten minutes ...
echo NSmGnprot : > "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/NSmGnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo NSmprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/NSmprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gnprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/Gnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gcprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/Gcprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo for the nucleic sequences...
echo this operation may take up to ten minutes ...
echo NSmGnnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmGnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo NSmnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmnucl_sdbl.aln-fasta -m MFP -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gnnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gcnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gcnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt


echo End of treatment for"${study}" :

#CCHFV
study="CCHFV"

echo Start of treatment for"${study}" : 

echo The position of the proteins of the "${study}" are being marked ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/mark_prot_CCHFV.py

echo The proteins are being splited for the "${study}" ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/split_protein_CCHFV.py

echo The NSm and Gn protein identifiers are currently being modified to avoid duplicate identifiers for the "${study}"...
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnprot.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnnucl.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta

mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl.aln-fasta

echo Non-informative regions are being removed and sequences that are too close will be deleted for the "${study}"...
./Phylogenetic/Script/trimal_CCHFV.sh

echo Identification of the best phylogenetic model is underway for the "${study}"

echo For the proteique sequences...
echo this operation may take up to ten minutes ...
echo NSmGnprot : > "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/NSmGnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo NSmprot : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/NSmprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo MLDprot : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/MLDprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo GP38prot : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/GP38prot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo Gnprot : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/Gnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo Gcprot : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
prottest -i "${path}"CCHFV/seg_M_phylo/Gcprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt


echo for the nucleic sequences...
echo this operation may take up to ten minutes ...
echo NSmGnnucl : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmGnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo NSmnucl : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo MLDnucl: >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/MLDnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo GP38nucl: >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/GP38nucl_sdbl.aln-fasta -m MFP -nt AUTO |grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo Gnnucl : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt

echo Gcnucl : >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gcnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"CCHFV/seg_M_phylo/Best_model_CCHFV.txt


echo End of treatment for"${study}" :

#Nairoviridae
study="Nairoviridae"

echo Start of treatment for"${study}" : 

echo The position of the proteins of the "${study}" are being marked ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/mark_prot_nairo.py

echo The position of the proteins of the "${study}" are being marked ...
/usr/bin/python3 /home/seannacharles/Phylogenetic/Script/split_protein_nairo.py

echo The NSm and Gn protein identifiers are currently being modified to avoid duplicate identifiers for the "${study}"...
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnprot.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnnucl.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot_align.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl_align.aln-fasta

echo Non-informative regions are being removed and sequences that are too close will be deleted for the "${study}"...
./Phylogenetic/Script/trimal_nairo.sh


echo Identification of the best phylogenetic model is underway for the "${study}"

echo For the proteique sequences...
echo this operation may take up to ten minutes ...
echo NSmGnprot : > "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/NSmGnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo NSmprot : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/NSmprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo MLDprot : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/MLDprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo GP38prot : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/GP38prot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo Gnprot : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/Gnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo Gcprot : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
prottest -i "${path}"Nairoviridae/seg_M_phylo/Gcprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo for the nucleic sequences...
echo this operation may take up to ten minutes ...
echo NSmGnnucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmGnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo NSmnucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo MLDnucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/MLDnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo GP38nucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/GP38nucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo Gnnucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo Gcnucl : >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gcnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}"Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt

echo End of treatment for"${study}" :

echo The phylogenetic trees are currently being created. This operation may take a few DAYS depending on the number of sequences. 
echo For this step you will need to install GNU Screen

# screen -dmS Phyml_scripts
# screen -S Phyml_scripts -X screen bash -c './Phylogenetic/Script/phyml_execution_CCHFV.sh'
# screen -S Phyml_scripts -X screen bash -c './Phylogenetic/Script/phyml_execution_Nairoviridae.sh'
# screen -S Phyml_scripts -X screen bash -c './Phylogenetic/Script/phyml_execution_Perybunyaviridae.sh'
# screen -ls
# while [ ! -f /home/seannacharles/Documents/scriptCCHFV_done.txt ] || [ ! -f /home/seannacharles/Documents/scriptNairoviridae_done.txt ] || [ ! -f /home/seannacharles/Documents/scriptPerybunaviridae_done.txt ]; do
#     sleep 1
# done

#screen -r Phyml_scripts
./Phylogenetic/Script/phyml_execution_CCHFV.sh && ./Phylogenetic/Script/phyml_execution_Nairoviridae.sh &&./Phylogenetic/Script/phyml_execution_Perybunyaviridae.sh &
#./Phylogenetic/Script/phyml_execution.sh

Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "It is finally over excution time : '${Hours}' : '${Minutes}' : '${Sec}' "