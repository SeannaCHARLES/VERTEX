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



for study,IO,protein,organismesAV in zip(names,INOUT,protINOUT,orgaWIWO) : 

    organismesAV=organismesAV.split(";")
    if IO !='' :
        organismes=[]
        for element in organismesAV : 
            organismes.append(Modif_ID(element[3:]))

        if IO.lower()=="in" :
         
            deleteIN(study,IO,protein,organismes,res_file,types_seq,path)
        else : 
            
            deleteOUT(study,IO,protein,organismes,res_file,types_seq,path) 