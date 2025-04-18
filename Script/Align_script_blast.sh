#!/bin/bash

#pwd
# /home/seannacharles

# makeblastdb    -in /home/seannacharles/Documents/Data/pdb_M.fasta \
#                 -title DB_M_seg \
#                 -dbtype nucl \
#                 -out /home/seannacharles/Documents/data_base/DB_M_seg \
#                 -parse_seqids

#BASE DE DONNÉE DÉJA CRÉE
#pdb_M.fasta contient toutes les séquences du segement M des virus suivant  : 


#Blast all de DB_M_seg alignement global
blastn  -query /home/seannacharles/Documents/Data/pdb_M.fasta \
       -db /home/seannacharles/Documents/data_base/DB_M_seg \
       -out /home/seannacharles/Documents/Align_herbevirus_ortho_Oropouch_global_all.txt \
       -outfmt 7 \
       -task blastn 

#Blast all de DB_M_seg alignement local
blastn  -query /home/seannacharles/Documents/Data/pdb_M.fasta \
       -db /home/seannacharles/Documents/data_base/DB_M_seg \
       -out /home/seannacharles/Documents/Align_herbevirus_ortho_Oropouch_local_all.txt \
       -outfmt 7  

#Blast Oropouch contre de DB_M_seg alignement global
blastn  -query /home/seannacharles/Documents/Data/Oropouch/M_KP052851.1_RNA.fna \
       -db /home/seannacharles/Documents/data_base/DB_M_seg \
       -out /home/seannacharles/Documents/Align_herbevirus_ortho_Oropouch_global.txt \
       -outfmt 7 \
       -task blastn 

#Blast Oropouch contre de DB_M_seg alignement local
blastn  -query /home/seannacharles/Documents/Data/Oropouch/M_KP052851.1_RNA.fna \
       -db /home/seannacharles/Documents/data_base/DB_M_seg \
       -out /home/seannacharles/Documents/Align_herbevirus_ortho_Oropouch_local.txt \
       -outfmt 7 \


       