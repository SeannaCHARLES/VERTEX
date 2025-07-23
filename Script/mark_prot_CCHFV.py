# Mark protein of CCHFV 
from Bio import SeqIO
from Bio.Seq import Seq

import subprocess
import sys 

import pandas as pd
from All_Functions import compt_ali_nucl,compt_alignement,mark_sequence


path =str(sys.argv[1])
study= str(sys.argv[2])
target=sys.argv[3]
stratseq=sys.argv[4]

posprot=sys.argv[5]
posprot=list(posprot.split(","))

all_prot_pos={}

df = pd.read_csv(f"{path}Pos_prot_studies.txt", delimiter="\t",header=0)
df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)
end_line=df.shape[0]







for line in range(1,end_line) : 
    if df.iloc[line,0] == study : 
        all_prot_pos[df.iloc[line,1]]=[df.iloc[line,2:]]

for record in SeqIO.parse(f"{path}{study}/prot_all_align_{study}.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_alignement(record.seq,pos_prot)

name=f"{path}{study}/prot_all_align_{study}_marke_prot.aln-fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()


for record in SeqIO.parse(f"{path}{study}/prot_all_align_{study}.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)
    



def compt_ali_nucl(sequence,debseq,pos_prot) :
    pos_align=[]
    compt_ali=0
    compt_nu=0
    compt_AA=0
    compt_codon=0
    codon=""
    deb=True

    for nu in sequence: 
        compt_ali+=1
        if nu!= "-" : 
            compt_nu+=1
            if compt_nu >=debseq: 
                compt_codon+=1
                codon+=nu
                if compt_codon==3 : 
                    compt_codon=0
                    compt_AA+=1
                    if compt_AA in pos_prot :
                            
                        if deb==True : 
                            deb=False
                            pos_align.append(compt_ali-2)
                            pos_prot.remove(compt_AA)
                        else : 
                            deb=True
                            pos_align.append(compt_ali)
                            pos_prot.remove(compt_AA) 
                        
                    codon=""
    return pos_align
  
##Function compt_alignement
# This function uses the position of the target organism's proteins (sequence proteic) to identify the position of the proteins in the alignment.

def compt_alignement (sequence,pos_prot) : 
    
    pos_align=[]
    compt_ali=0
    compt_AA=0
    mark_prot=""
    deb=True
    for aa in sequence : 
        if aa=='-' : 
            compt_ali+=1
        else : 
            compt_ali+=1
            compt_AA+=1

        if compt_AA in pos_prot :
            if deb==True : 
                mark_prot+="|"
                mark_prot+=aa
                deb=False
                pos_align.append(compt_ali)
                pos_prot.remove(compt_AA)
            else : 
                mark_prot+=aa
                mark_prot+="|"
                pos_align.append(compt_ali)
                deb=True
                pos_prot.remove(compt_AA)
        else : 
            mark_prot+=aa
    return pos_align

##Function mark_sequence
#This function marks the start and end of proteins in protein and nucleic acid sequences. 

def mark_sequence(sequence,pos_align,name,description) : 
    compt_ali=0
    mark_prot=""
    deb=True
    for aa in sequence : 
        compt_ali+=1
        if compt_ali in pos_align :
            if deb==True : 
                mark_prot+="|"
                mark_prot+=aa
                deb=False
             
            else : 
                mark_prot+=aa
                mark_prot+="|"
                deb=True
        else : 
            mark_prot+=aa
    
    file = open(name,"a")
    text=str (">"+description +'\n'+mark_prot+'\n')
    file.write(text)
    file.close()



















#prot

#nucl
for record in SeqIO.parse(f"{path}{study}/nucl_all_align_{study}.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_ali_nucl(record.seq,stratseq,pos_prot)
              
name=f"{path}{study}/nucl_all_align_{study}_marke.aln-fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}{study}/nucl_all_align_{study}.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)

