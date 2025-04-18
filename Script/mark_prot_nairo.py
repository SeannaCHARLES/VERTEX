from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
from All_Function import compt_ali_nucl,compt_alignement,mark_sequence
path ="/home/seannacharles/Documents/Nairoviridae/"
target = "haemorrha"
print('mark prot nairoviridae')
pos_align=[]

debseq=93
# def compt_ali_nucl(sequence) :
#     pos_align=[]
#     compt_ali=0
#     compt_nu=0
#     compt_prot=0
#     compt_codon=0
#     codon=""
#     deb=True
#     for nu in record.seq : 
#         compt_ali+=1
#         if nu!= "-" : 
#             compt_nu+=1
#             if compt_nu >=93: 
#                 compt_codon+=1
#                 codon+=nu
#                 if compt_codon==3 : 
#                     compt_codon=0
#                     compt_prot+=1
#                     if compt_prot in pos_prot :
                            
#                         if deb==True : 
#                             deb=False
#                             pos_align.append(compt_ali-2)
#                             pos_prot.remove(compt_prot)
#                         else : 
#                             deb=True
#                             pos_align.append(compt_ali)
#                             pos_prot.remove(compt_prot) 
                        
#                     codon=""
#     return pos_align

# def compt_alignement (sequence) : 
#     global pos_prot
#     pos_align=[]
#     compt_ali=0
#     compt_prot=0
#     mark_prot=""
#     deb=True
#     for aa in sequence : 
#             if aa=='-' : 
#                 compt_ali+=1
#             else : 
#                 compt_ali+=1
#                 compt_prot+=1

#             if compt_prot in pos_prot :
#                 if deb==True : 
#                     mark_prot+="|"
#                     mark_prot+=aa
#                     deb=False
#                     pos_align.append(compt_ali)
#                     pos_prot.remove(compt_prot)
#                 else : 
#                     mark_prot+=aa
#                     mark_prot+="|"
#                     pos_align.append(compt_ali)
#                     deb=True
#                     pos_prot.remove(compt_prot)
#             else : 
#                 mark_prot+=aa
#     return pos_align

# def mark_sequence(sequence,pos_align,name) : 
#     compt_ali=0
#     mark_prot=""
#     deb=True
#     for aa in sequence : 
#         compt_ali+=1
#         if compt_ali in pos_align :
#             if deb==True : 
#                 mark_prot+="|"
#                 mark_prot+=aa
#                 deb=False
             
#             else : 
#                 mark_prot+=aa
#                 mark_prot+="|"
#                 deb=True
#         else : 
#             mark_prot+=aa
#     file = open(name,"a")
#     text=str (">"+record.description +'\n'+mark_prot+'\n')
#     file.write(text)
#     file.close()

# prot
pos_prot=[23,248,249,520,521,843,844,996,997,1684]

for record in SeqIO.parse(f"{path}M_prot_all_align_Nairoviridae.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_alignement(record.seq,pos_prot)



name=f"{path}M_prot_all_align_marke_prot.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}M_prot_all_align_Nairoviridae.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


# nucl
pos_prot=[23,248,249,520,521,843,844,996,997,1684]

for record in SeqIO.parse(f"{path}M_nucl_all_align_Nairoviridae.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_ali_nucl(record.seq,debseq,pos_prot)



name=f"{path}M_nucl_all_align_marke.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}M_nucl_all_align_Nairoviridae.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


#subprocess.run("/usr/bin/python3 /home/seannacharles/Documents/Script/split_protein_nairo.py",shell=True)