# Script qui sépare les différentes proteines du segement M
from Bio import SeqIO
import re
# IL FAUT SUPPRIMER LES FICHIERS CRÉE PAR LE SCRIPT AVANT DE LE RELANCER
AA_L=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","-","|"]
Protname=["MLD_M","GP38_M","Gn_M","NSm_M","Gc_M"]
vi=["M","L","S"]
Data={}

for segment in vi : 
    for record in SeqIO.parse(f"/home/seannacharles/Documents/Perybunyaviridae/{segment}_prot_all_align_X.fst","fasta") : 
        if record.id not in Data :
            seq=re.findall(r'\|(.*?)\|',str (record.seq))
            for index,prot in enumerate(seq) : 
                for AA in prot:
                    if AA not in AA_L : 
                        prot_mod=prot.replace(str(AA),"N")
                        seq[index]=prot_mod
                    
            Data[record.id]=seq
    

    for ID,seq in Data.items() : 
        for index,prot in enumerate(seq) : 
            name="/home/seannacharles/Documents/Perybunyaviridae/seg_M_phylo"+str(Protname[index])+".fasta"
            file = open(name,"a")
            text=str (">"+ID+'\n'+prot+'\n')
            file.write(text)
            file.close()



