from Bio import SeqIO
import re
import subprocess
from Bio.Seq import Seq
from Bio import AlignIO
from All_Function import Extract_data, write_file
path = "/home/seannacharles/Documents/Nairoviridae/"
study="Nairoviridae"

Protname=["MLD","GP38","Gn","NSm","Gc"]


# def Extract_data(file,sp) : 
#     if "nucl" in sp  :
#         AA_L=["a","t","c","g","-"]
#     else :
#         AA_L=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","-","|"]
#     Data={}

#     for record in SeqIO.parse(f"{path}{file}","fasta") :
#         if record.id not in Data :
#             seq=re.findall(r'\|(.*?)\|',str (record.seq))
#             for index,prot in enumerate(seq) : 
#                 for AA in prot:
#                     if AA not in AA_L : 
#                         prot_mod=prot.replace(str(AA),"n")
#                         seq[index]=prot_mod
                    
#             Data[record.description]=seq
#     write_file_prot(Data,sp)
    
# def write_file_prot(Data,sp) : 
#     global Protname
 
#     for ID,seq in Data.items() : 
#         for index,prot in enumerate(seq) : 
#             name=f"{path}seg_M_phylo/"+str(Protname[index])+sp+".fasta"
#             file = open(name,"a")
#             text=str (">"+ID+'\n'+Seq(prot.upper())+'\n')
#             file.write(text)
#             file.close()
            
            

print('split prot Nairoviridae')



Extract_data("M_prot_all_align_marke_prot.fasta","prot",Protname,study)

Extract_data("M_nucl_all_align_marke.fasta","nucl",Protname,study)





# subprocess.run(f"cat {path}seg_M_phylo/NSmprot.fasta {path}seg_M_phylo/Gnprot.fasta > {path}seg_M_phylo/NSmGnprot.fasta",shell=True)
# subprocess.run(f"cat {path}seg_M_phylo/NSmnucl.fasta {path}seg_M_phylo/Gnnucl.fasta > {path}seg_M_phylo/NSmGnnucl.fasta",shell=True)

subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnprot.fasta > {path}seg_M_phylo/Gnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcprot.fasta > {path}seg_M_phylo/Gcprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmprot.fasta > {path}seg_M_phylo/NSmprot_align.aln-fasta",shell=True)
#subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmGnprot.fasta > {path}seg_M_phylo/NSmGnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38prot.fasta > {path}seg_M_phylo/GP38prot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDprot.fasta > {path}seg_M_phylo/MLDprot_align.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnnucl.fasta > {path}seg_M_phylo/Gnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcnucl.fasta > {path}seg_M_phylo/Gcnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmnucl.fasta > {path}seg_M_phylo/NSmnucl_align.aln-fasta",shell=True)
#subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmGnnucl.fasta > {path}seg_M_phylo/NSmGnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38nucl.fasta > {path}seg_M_phylo/GP38nucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDnucl.fasta > {path}seg_M_phylo/MLDnucl_align.aln-fasta",shell=True)

#subprocess.run("./Documents/Script/trimal_nairo.sh")