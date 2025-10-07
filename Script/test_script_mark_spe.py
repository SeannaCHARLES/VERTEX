
##########################Description##########################
""" This script seperate the different protein according to their 
real position or the organism of interest position.

"""
##########################################################

##########################Import##########################
from Bio import SeqIO
from Bio.Seq import Seq

import subprocess
import sys 
import pandas as pd
from All_Functions import nucleic_sequences,proteic_sequence,final_pos,col_pos
##########################################################

###########################Input##########################
path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
position=list(sys.argv[3].split(" "))
names=list(sys.argv[4].split(" "))
positionspeOG=list(sys.argv[5].split(" "))
start=list(sys.argv[6].split(" "))
startspeOG=list(sys.argv[7].split(" "))
prot=list(sys.argv[8].split(" "))
protspe_OG=list(sys.argv[9].split(" "))
resfile=str(sys.argv[10])
studiesspeOG=list(sys.argv[11].split(" "))
NUIDspeOG=list(sys.argv[12].split(" "))

onlyOG= [x for x in names if x not in studies]
Pos_prot={}
proteins={}
study_prot={}
##########################################################

###########################Script##########################
for study,start,position in zip(studies,start,position) :
    position=position.split(";")
    position=[int(x) for x in position]
    info =  position+[int(start)]
    proteins[study]=info

for study,start,position in zip(onlyOG,startspeOG,positionspeOG) :
    position=position.split(";")
    position=[int(x) for x in position]
    info =  position+[int(start)]
    proteins[study]=info


for item,study in zip(protspe_OG,onlyOG) :
        study_prot[study]=item.split(";")

for item,study in zip(prot,studies) :
        study_prot[study]=item.split(";")



for study in studies :
    subprocess.run(["rm",f"{path}{study}{resfile}*" ],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,shell=True)
    final_pos_prot,name=final_pos(study,prot,study_prot,proteins)
    sp=study
    file = open(f"{name}","a")

    for key,positions in final_pos_prot.items() : 
        file.write(f"{study}\t{key}\t")

        for pos in positions : 
            file.write(f"{pos}\t")
        file.write(f"\n")

    file.close()

    for record in SeqIO.parse(f"{path}{study}/prot_sdbl_{study}.fasta","fasta") :
        proteic_sequence(study,record.id,record.seq,final_pos_prot[record.id],study_prot[study],sp,resfile)


    for record in SeqIO.parse(f"{path}{study}/nucl_sdbl_{study}.fasta","fasta") :
        stratseq=proteins[sp][-1]
        nucleic_sequences(study,record.id,record.seq,stratseq,final_pos_prot[record.id],study_prot[study],sp,resfile)

for study,sp,ID in zip(studiesspeOG,onlyOG,NUIDspeOG) :

    final_pos_prot,name=final_pos(study,prot,study_prot,proteins)
    file = open(f"{name}","a")

    for key,positions in final_pos_prot.items() : 
        file.write(f"{study}\t{key}\t")

        for pos in positions : 
            file.write(f"{pos}\t")
        file.write(f"\n")

    file.close()

    for record in SeqIO.parse(f"{path}{study}/{ID}_prot.fasta","fasta") :
        proteic_sequence(study,record.id,record.seq,final_pos_prot[record.id],study_prot[study],sp,resfile)


    for record in SeqIO.parse(f"{path}{study}/{ID}_nucl.fasta","fasta") :
        stratseq=proteins[sp][-1]
        nucleic_sequences(study,record.id,record.seq,stratseq,proteins[sp][:-1],study_prot[study],sp,resfile)
##########################################################
