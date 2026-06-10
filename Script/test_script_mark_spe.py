
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

from All_Functions import nucleic_sequences,proteic_sequence,final_pos,col_pos,final_pos_IO,nucleic_sequences_OG
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
orgaWIWO=list(sys.argv[13].replace(" ", "\n").splitlines())
INOUT=list(sys.argv[14].replace(" ", "\n").splitlines())
protINOUT=list(sys.argv[15].replace(" ", "\n").splitlines())
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

concern_study={}
for name_,orga,IO,proteinIO in zip(names,orgaWIWO,INOUT,protINOUT) : 
    if orga!="" : 
        orga=orga.split(";")
        conc_orga=[element.replace(".", "").lower() for element in orga]
        
        concern_study[name_]=[conc_orga,IO,proteinIO]
for study in studies : 
#for name_all,organismesAV,IO,proteinIO in zip(names,orgaWIWO,INOUT,protINOUT) :
    #subprocess.run(["rm",f"{path}{study}{resfile}/*" ],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,shell=True)
   
    # if study in concern_study : 
    #     final_pos_prot,name=final_pos_IO(study,prot,study_prot,proteins,concern_study)
    #     #final_pos_prot,name=final_pos(study,prot,study_prot,proteins)
    # else : 
    #     final_pos_prot,name=final_pos(study,prot,study_prot,proteins)

    sp=study
    # file = open(f"{name}","a")

    # for key,positions in final_pos_prot.items() : 
    #     file.write(f"{study}\t{key}\t")

    #     for pos in positions : 
    #         file.write(f"{pos}\t")
    #     file.write(f"\n")

    # file.close()
    
    for record in SeqIO.parse(f"{path}{study}/prot_sdbl_{study}.fasta","fasta") :
        proteic_sequence(study,record.id,record.seq,study_prot[study],sp,resfile,path)


    for record in SeqIO.parse(f"{path}{study}/nucl_sdbl_{study}.fasta","fasta") :
        stratseq=proteins[sp][-1]
        nucleic_sequences(study,record.id,record.seq,stratseq,study_prot[study],sp,resfile,path)

for study,sp,OG_list,name,organismesAV,IO,proteinIO in zip(studiesspeOG,onlyOG,NUIDspeOG,names,orgaWIWO,INOUT,protINOUT) :

    # if study in concern_study : 
    #     #final_pos_prot,name=final_pos_(study,prot,study_prot,proteins,concern_study)
    #     final_pos_prot,name=final_pos(study,prot,study_prot,proteins)
    # else : 
    #     final_pos_prot,name=final_pos(study,prot,study_prot,proteins)
    # file = open(f"{name}","a")

    # for key,positions in final_pos_prot.items() : 
    #     file.write(f"{study}\t{key}\t")

    #     for pos in positions : 
    #         file.write(f"{pos}\t")
    #     file.write(f"\n")

    # file.close()
    OG=OG_list.split(";")
    for ID in OG : 
        for record in SeqIO.parse(f"{path}{study}/{ID}_prot.fasta","fasta") :
        
            proteic_sequence(study,record.id,record.seq,study_prot[study],sp,resfile,path)


        for record in SeqIO.parse(f"{path}{study}/{ID}_nucl.fasta","fasta") :
            stratseq=proteins[sp][-1]
        
            nucleic_sequences_OG(study,record.id,record.seq,stratseq,proteins[sp][:-1],study_prot[study],sp,resfile,path)
##########################################################
