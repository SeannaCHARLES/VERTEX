from Bio import SeqIO
from Bio.Seq import Seq

import subprocess

from All_Function import compt_ali_nucl,compt_alignement,mark_sequence

path ="/home/seannacharles/Documents/Perybunyaviridae/"

debseq=57
target ="bunyamwer"

#All peri prot
pos_prot=[22,302,303,476,509,1372]

for record in SeqIO.parse(f"{path}M_prot_all_align_.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_alignement(record.seq,pos_prot)


name=f"{path}M_prot_all_align_marke_prot.aln-fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}M_prot_all_align_.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


#All peri nucl
pos_prot=[22,302,303,476,509,1372]

for record in SeqIO.parse(f"{path}M_nucl_all_align_.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_ali_nucl(record.seq,debseq,pos_prot)



name=f"{path}M_nucl_all_align_marke.aln-fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}M_nucl_all_align_.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)



pos_prot=[22,302,303,476,509,1372]

#hedwig prot
for record in SeqIO.parse(f"{path}M_prot_all_align_hedwig.aln-fasta","fasta"):
    if target in record.description:
        pos_align=compt_alignement(record.seq,pos_prot)

name=f"{path}M_prot_all_align_marke_prot_hedwig.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(name,"w")
file.close()

for record in SeqIO.parse(f"{path}M_prot_all_align_hedwig.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)
pos_prot=[22,302,303,476,509,1372]

#hedwig nucl
for record in SeqIO.parse(f"{path}M_nucl_all_align_hedwig.aln-fasta","fasta"):
    if target in record.description:
        pos_align=compt_ali_nucl(record.seq,debseq,pos_prot)
        
name=f"{path}M_nucl_all_align_marke_hedwig.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(name,"w")
file.close()

for record in SeqIO.parse(f"{path}M_nucl_all_align_hedwig.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


#khurdun prot
pos_prot=[23,292,526,944]


for record in SeqIO.parse(f"{path}LEIV-Ast01-5_YP_010840266.1_khurdunyp.fasta","fasta") :
    pos_align=compt_alignement(record.seq,pos_prot)
    
name=f"{path}M_prot_marke_prot_khurdun.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(name,"w")
file.close()

debseq=34
for record in SeqIO.parse(f"{path}LEIV-Ast01-5_YP_010840266.1_khurdunyp.fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)
pos_prot=[23,292,526,944]   

#khurdun nucl
for record in SeqIO.parse(f"{path}LEIV-Ast01-5_NC_078560.1_khurdunyp.fasta","fasta") :
    pos_align=compt_ali_nucl(record.seq,debseq,pos_prot)
    
name=f"{path}M_nucl_marke_khurdun.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(name,"w")
file.close()


for record in SeqIO.parse(f"{path}LEIV-Ast01-5_NC_078560.1_khurdunyp.fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)
