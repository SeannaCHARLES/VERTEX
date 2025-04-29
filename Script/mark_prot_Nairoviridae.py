from Bio import SeqIO
from Bio.Seq import Seq

import subprocess

from All_Functions import compt_ali_nucl,compt_alignement,mark_sequence

path ="/home/seannacharles/Documents/"
study="Nairoviridae"

target = "haemorrha"

pos_align=[]
stratseq=93

# prot
pos_prot=[23,248,249,520,521,843,844,996,997,1684]

for record in SeqIO.parse(f"{path}{study}/M_prot_all_align_{study}.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_alignement(record.seq,pos_prot)



name=f"{path}{study}/M_prot_all_align_marke_prot.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}{study}/M_prot_all_align_{study}.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


# nucl
pos_prot=[23,248,249,520,521,843,844,996,997,1684]

for record in SeqIO.parse(f"{path}{study}/M_nucl_all_align_{study}.aln-fasta","fasta") :
    if target in record.description:
        pos_align=compt_ali_nucl(record.seq,stratseq,pos_prot)



name=f"{path}{study}/M_nucl_all_align_marke.fasta"
subprocess.run(f"rm {name}",shell=True)
file = open(f"{name}","w")
file.close()

for record in SeqIO.parse(f"{path}{study}/M_nucl_all_align_{study}.aln-fasta","fasta") :
    mark_sequence(record.seq,pos_align,name,record.description)


