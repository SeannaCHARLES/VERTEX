
import subprocess
import sys 
import numpy as np

from All_Functions import first_call,concatenates,extract_id
from itertools import chain

path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
studiesspeOG=list(sys.argv[3].split(" "))

names=list(sys.argv[4].split(" "))
NUIDspeOG=list(sys.argv[5].split("''"))[0].split("  ")

email=sys.argv[6]

list_prot = [x.split(';') for x in sys.argv[7].split()]
list_prot_speOG = [x.split(';') for x in sys.argv[8].split()]

positionspeOG=list(sys.argv[8].split(" "))
OG_prot=list(sys.argv[9].split(" "))

orgaWIWO=list(sys.argv[10].replace(" ", "\n").splitlines())
INOUT=list(sys.argv[11].replace(" ", "\n").splitlines())
protINOUT=list(sys.argv[12].replace(" ", "\n").splitlines())

Allstudies=list(sys.argv[13].split(" "))
ALLNUID=list(sys.argv[14].split(" "))

namespeprot=list(sys.argv[15].split(" "))
namespeprotOG=list(sys.argv[16].split(" "))





for pos,element in enumerate(namespeprot) : 
    namespeprot[pos]=element.split(";")
for pos,element in enumerate(namespeprotOG) : 
    namespeprotOG[pos]=element.split(";")

chain_prot=""

dicy={}
vv_prot={}

for all_s, IO,prots,virus,list_ID in zip(Allstudies,INOUT,protINOUT,orgaWIWO,ALLNUID) :
   
    list_ID=list_ID.split(";")
 
    
    if IO!='' :
        if all_s not in dicy : 
            dicy[all_s]={}
       
        if IO=="OUT" : 
            for v in virus.split(";") : 
                if v not in dicy[all_s] :
                    dicy[all_s][v]=[prots] 
                else : 
                    dicy[all_s][v].append([prots])

        elif IO=="IN" : 
            for  ID in list_ID :
                if ID not in virus.split(";") :
                    if ID not in dicy[all_s] :
                        dicy[all_s][ID]=[prots]
                    elif prots not in  dicy[all_s][ID] : 
                        dicy[all_s][ID].append(prots)
        elif IO=='none' :
            if all_s not in dicy : 
                dicy[all_s]={}

if not dicy : 
    for all_s in Allstudies : 
        dicy[all_s]=["-","-","-"]





for study,allprot,spe_prot in zip(studies,list_prot,namespeprot) :
   
    chain_prot=""
    print(f"For the study : {study}")
    subprocess.run(f"rm {path}{study}/nucl/* {path}{study}/prot/* 2>/dev/null || true",shell=True)
    
    prot=list(dict.fromkeys(allprot))
    vv_prot[study]=prot
    
    for p in prot : 

        chain_prot+=f"Start{p}\tEnd{p}\t"

    subprocess.run(f"echo 'Study\tVirus name\tID itol tree\tTaxonomy\tInformations\tProtein accesion number\tSize Proteic sequence\t{chain_prot}Nucleic accession number\tSize Nucleic sequence\tStrand\tStrain\tGeographical region\tHost or vector\tReferences' >{path}{study}/metadata_{study}.txt",shell=True)
    
    subprocess.run(f"echo 'Study\tID itol tree\t{chain_prot}' >{path}{study}/pos_prot_{study}.txt",shell=True)

    with open(f"{path}{study}/sequences_{study}.gb","r") as Datas:
        sp=''
        first_call(path,study,spe_prot)
        
        extract_id(path,Datas,study,sp,prot,email,vv_prot[study],dicy,spe_prot)
        concatenates(path,study)


    print(f"Mafft alignement for the file {path}{study}/prot_all_{study}.fasta save in the output file : {path}{study}/prot_all_align_{study}.aln-fasta IN PROGRESS...")
    subprocess.run(f"echo Mafft for the outfile {path}{study}/prot_all_align_{study}.aln-fasta >> {path}{study}/Maffts.log && mafft --auto {path}{study}/prot_all_{study}.fasta > {path}{study}/prot_all_align_{study}.aln-fasta 2>> {path}{study}/Maffts.log",shell=True)

    print(f"Mafft alignement for the file {path}{study}/nucl_all_{study}.fasta save in the output file : {path}{study}/nucl_all_align_{study}.aln-fasta IN PROGRESS...")
    subprocess.run(f"echo Mafft for the outfile {path}{study}/nucl_all_align_{study}.aln-fasta >> {path}{study}/Maffts.log && mafft --auto {path}{study}/nucl_all_{study}.fasta > {path}{study}/nucl_all_align_{study}.aln-fasta 2>> {path}{study}/Maffts.log",shell=True)

for study,OG_list,prot,spe_prot in zip(studiesspeOG, NUIDspeOG,list_prot_speOG,namespeprotOG):
   
    OG=OG_list.split(" ")
    for ID in OG : 
        if ID!="" : 
            with open(f"{path}{study}/{ID}.gb","r") as Datas:
                sp='OGspe'
                extract_id(path,Datas,study,sp,prot,email,vv_prot[study],speci_order,dicy,spe_prot)
