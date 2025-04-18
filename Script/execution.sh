#!/bin/bash
SECONDS=0
path="/home/seannacharles/Documents/"

#./Documents/Script/download_sequences.sh
/usr/bin/python3 "${path}"Script/genbank_extract.py

Perybunyaviridae
echo Start Perybunyaviridae
study="Perybunyaviridae"
/usr/bin/python3 "${path}"Script/mark_prot.py
/usr/bin/python3 "${path}"Script/split_protein_peri.py

{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprotSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnprotSH.fasta ; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnuclSH.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclKhurdun.fasta;sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/GnnuclSH.fasta ;} > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot_align.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl_align.aln-fasta
./Documents/Script/trimal_peri.sh

echo NSmGnprot : > "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/NSmGnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo NSmprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/NSmprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gnprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/Gnprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gcprot : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
prottest -i "${path}""${study}"/seg_M_phylo/Gcprot_sdbl.aln-fasta -all  -DT -F |grep 'Best model according to DT:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt


echo NSmGnnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmGnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo NSmnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/NSmnucl_sdbl.aln-fasta -m MFP -nt AUTO |grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gnnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gnnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt

echo Gcnucl : >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt
iqtree2 -s "${path}""${study}"/seg_M_phylo/Gcnucl_sdbl.aln-fasta -m MFP -nt AUTO|grep 'Best-fit model:' >> "${path}""${study}"/seg_M_phylo/Best_model_Perybunyaviridae.txt


echo End Perybunyaviridae

#CCHFV
echo Start CCHFV
study="CCHFV"
/usr/bin/python3 "${path}"Script/mark_prot_CCHFV.py
/usr/bin/python3 "${path}"Script/split_protein_CCHFV.py
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnprot.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnnucl.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl.aln-fasta

./Documents/Script/trimal_CCHFV.sh

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


echo End CCHFV

#Nairoviridae
echo Start Nairoviridae
study="Nairoviridae"
/usr/bin/python3 "${path}"Script/mark_prot_nairo.py
/usr/bin/python3 "${path}"Script/split_protein_nairo.py
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmprot.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnprot.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta
{ sed 's/^>/>N/' "${path}""${study}"/seg_M_phylo/NSmnucl.fasta; sed 's/^>/>G/' "${path}""${study}"/seg_M_phylo/Gnnucl.fasta; } > "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnprot.fasta > "${path}""${study}"/seg_M_phylo/NSmGnprot_align.aln-fasta
mafft --auto "${path}""${study}"/seg_M_phylo/NSmGnnucl.fasta > "${path}""${study}"/seg_M_phylo/NSmGnnucl_align.aln-fasta

./Documents/Script/trimal_nairo.sh

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

echo End Nairoviridae

Hours=$(("${SECONDS}" / 3600))
Minutes=$((("${SECONDS}" % 3600) / 60))
Sec=$(("${SECONDS}" % 60))

echo "execution time : '${Hours}' : '${Minutes}' : '${Sec}' "