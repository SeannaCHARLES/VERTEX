# translate prot seq in nu for selective pressure 
from Bio import SeqIO
from Bio.Seq import *

pos_prot=[23,248,249,520,521,843,844,996,997,1684]
pos_nu=[]
name="/home/seannacharles/Documents/CCHFV/M_nucl_all_align_marke_nucl.aln-fasta" 


target="NC_005300.2"


for record in SeqIO.parse("/home/seannacharles/Documents/CCHFV/M_nucl_all_align.aln-fasta","fasta") :
    if target in record.name:
        compt_ali=0
        compt_nu=0
        compt_prot=0
        compt_cod=1
        deb=True
        for nu in record.seq : 
            compt_ali+=1
            if nu!='-' : 
                
                if compt_nu >=92 : 
                    
                    if compt_cod==3 : 
                        compt_prot+=1


                        if compt_prot in pos_prot :
                            if deb==True : 
                                deb=False
                                pos_nu.append(compt_ali-2)
                                pos_prot.remove(compt_prot)
                            else : 
                                
                                pos_nu.append(compt_ali)
                                deb=True
                                pos_prot.remove(compt_prot)
                        
                            
                        
                        compt_cod=1
                    else : 
                        compt_cod+=1
                compt_nu+=1
            
                
                

for record in SeqIO.parse("/home/seannacharles/Documents/CCHFV/M_nucl_all_align.aln-fasta","fasta") :
    compt_ali=0
    mark_prot=""
    deb=True
    for aa in record.seq : 
        compt_ali+=1
        if compt_ali in pos_nu :
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
    text=str (">"+record.description +'\n'+Seq(mark_prot)+'\n')
    file.write(text)
file.close()

pos_prot=[23,248,249,520,521,843,844,996,997,1684]

pos_align=[]

for record in SeqIO.parse("/home/seannacharles/Documents/CCHFV/M_prot_all_align.aln-fasta","fasta") :
    if target in record.name:
        compt_ali=0
        compt_prot=0
        mark_prot=""
        deb=True
        for aa in record.seq : 
            if aa=='-' : 
                compt_ali+=1
            else : 
                compt_ali+=1
                compt_prot+=1

            if compt_prot in pos_prot :
                if deb==True : 
                    mark_prot+="|"
                    mark_prot+=aa
                    deb=False
                    pos_align.append(compt_ali)
                    pos_prot.remove(compt_prot)
                else : 
                    mark_prot+=aa
                    mark_prot+="|"
                    pos_align.append(compt_ali)
                    deb=True
                    pos_prot.remove(compt_prot)
            else : 
                mark_prot+=aa



name="/home/seannacharles/Documents/CCHFV/M_prot_all_align_marke_prot.aln-fasta"

for record in SeqIO.parse("/home/seannacharles/Documents/CCHFV/M_prot_all_align.aln-fasta","fasta") :
    compt_ali=0
    mark_prot=""
    deb=True
    for aa in record.seq : 
        compt_ali+=1
        if compt_ali in pos_nu :
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
    text=str (">"+record.description +'\n'+Seq(mark_prot)+'\n')
    file.write(text)
file.close()