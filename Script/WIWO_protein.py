import sys
import subprocess

from Bio import SeqIO
from All_Functions import Modif_ID,deleteIN,deleteOUT

path = sys.argv[1]
names=list(sys.argv[2].split(" "))
orgaWIWO=list(sys.argv[3].replace(" ", "\n").splitlines())
INOUT=list(sys.argv[4].replace(" ", "\n").splitlines())
protINOUT=list(sys.argv[5].replace(" ", "\n").splitlines())
res_file=str(sys.argv[6])
types_seq=list(sys.argv[7].split(";"))


def deleteIN(study,IO,protein,organismes,res_file,types_seq) : 

    for type_seq in types_seq : 
        subprocess.run(f" cp {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta {path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta" ,shell=True)
        subprocess.run(f"rm {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta",shell=True)
        filename=f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"

        for record in SeqIO.parse(f"{path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta","fasta") :
            print(f"organismes : {org[:6]} \n Record : {record.id}")
            if any(org[:6] in record.id for org in organismes):
                print ("write")
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()

def deleteOUT(study,IO,protein,organismes,res_file,types_seq) : 

    for type_seq in types_seq : 
        subprocess.run(f"cp {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta {path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta",shell=True)
        subprocess.run(f"rm {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta",shell=True)
        filename=f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"

        for record in SeqIO.parse(f"{path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta","fasta") :
            if not any(org[:6] in record.id for org in organismes): 
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()


for study,IO,protein,organismesAV in zip(names,INOUT,protINOUT,orgaWIWO) : 

    organismesAV=organismesAV.split(";")
    if IO !='' :
        organismes=[]
        for element in organismesAV : 
            organismes.append(Modif_ID(element[3:]))

        if IO.lower()=="in" :
            print(f"for the study {study} the protein {protein} is {IO} the following organismes only {organismes}")
            deleteIN(study,IO,protein,organismes,res_file,types_seq)
        else : 
            print(f"for the study {study} the protein {protein} is {IO} the following organismes only {organismes}")
            deleteOUT(study,IO,protein,organismes,res_file,types_seq) 