
import sys
import subprocess
from All_Functions import select_doublon,delete_dbl,col_name
from Bio import SeqIO
path =str(sys.argv[1])
study=str(sys.argv[2])
target=str(sys.argv[3])
ths=float(sys.argv[4])
Prot_list=sys.argv[5]
Prot_list=list(Prot_list.split(","))
res_file=str(sys.argv[6])
OG_list=sys.argv[7]
OG_list=list(OG_list.split(","))
# study="Nairoviridae"
# target="haemorrha"
# ths=0.005
# Prot_list=["MLD","GP38","Gn","NSm","Gc"]
# path="/home/seannacharles/"
# res_file="/res_phylo/"

final_dbl=[]
old_col_names=[]

for protein in Prot_list : 
    Doublon,name_colonnes=select_doublon(protein,study,target,ths,OG_list)
    old_col_names=col_name(old_col_names,name_colonnes,protein)

    if len(final_dbl)==0 : 
        final_dbl=Doublon
    else : 
        final_dbl += [e for e in Doublon if e not in final_dbl]

for protein in Prot_list : 
        delete_dbl(final_dbl,study,protein,old_col_names)

doublon_file=f"{path}{study}/Doublon_file_{study}.txt"
dbl_file=open(doublon_file,'a')

for db in final_dbl : 
    dbl_file.write(f"{db}\n")
dbl_file.close()


NG_old_col_names=[]
NG_final_dbl=[]

for name in old_col_names:
    NG_old_col_names.append(str("N"+name))
    NG_old_col_names.append(str("G"+name))

NG_old_col_names.append(str("Gkhurdunyp"))#Khurdun is not present in the sequence file for the other proteins because it will only be used for NSmGn as an outgroup.

for name in final_dbl:
    NG_final_dbl.append(str("N"+name))
    NG_final_dbl.append(str("G"+name))

#Only the orgarnism in the NSm_listeNairo have NSm the other NSm sequences must be deleted
delete_dbl(NG_final_dbl,study,"NSmGn",NG_old_col_names)
NSm_listeNairo=["dugbeense","haemorrha","kupeyp010","ganjamamt"]
Concern_prot=["NSm","NSmGn"]
NUPR=["nucl","prot"]
if study=="Nairoviridae" : 
    for types in NUPR : 
        subprocess.run(f"cp {path}{study}/{res_file}NSm{types}_sdbl.fasta {path}{study}/{res_file}NSm{types}_sdbl_RN.fasta ",shell=True)
        subprocess.run(f"rm {path}{study}/{res_file}NSm{types}_sdbl.fasta",shell=True)
        filename=f"{path}{study}/{res_file}NSm{types}_sdbl.fasta"
        for record in SeqIO.parse(f"{path}{study}/{res_file}NSm{types}_sdbl_RN.fasta","fasta") :
            if record.id in NSm_listeNairo : 
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()

        subprocess.run(f"cp {path}{study}/{res_file}NSmGn{types}_sdbl.fasta {path}{study}/{res_file}NSmGn{types}_sdbl_RN.fasta ",shell=True)
        subprocess.run(f"rm {path}{study}/{res_file}NSmGn{types}_sdbl.fasta",shell=True)
        filename=f"{path}{study}/{res_file}NSmGn{types}_sdbl.fasta"
        for record in SeqIO.parse(f"{path}{study}/{res_file}NSmGn{types}_sdbl_RN.fasta","fasta") :
            if record.id.startswith("N") :
                if record.id[1:] in NSm_listeNairo: 
                    file = open(filename,"a")
                    file.write(f">{record.id}\n{record.seq}\n")       
                    file.close()
            else : 
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()