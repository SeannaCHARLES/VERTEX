from Bio import SeqIO
import re
import subprocess
from Bio.Seq import Seq
from Bio import AlignIO
from All_Functions import Extract_data, write_file
path = "/home/seannacharles/Documents/Peribunyaviridae/"
study="Peribunyaviridae"


 

print('split prot Peribunyaviridae')


Protname=["Gn","NSm","Gc"]
Extract_data("M_prot_all_align_marke_prot.aln-fasta","protSH",Protname,study)
Extract_data("M_nucl_all_align_marke.aln-fasta","nuclSH",Protname,study)

Extract_data("M_prot_all_align_marke_prot_hedwig.fasta","prot",Protname,study)
Extract_data("M_nucl_all_align_marke_hedwig.fasta","nucl",Protname,study)

Protname=["Gn","Gc"]
Extract_data("M_prot_marke_prot_khurdun.fasta","protKhurdun",Protname,study)
Extract_data("M_nucl_marke_khurdun.fasta","nuclKhurdun",Protname,study)






subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/Gnprot.fasta > {path}{study}seg_M_phylo/Gnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/Gcprot.fasta > {path}{study}seg_M_phylo/Gcprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/NSmprot.fasta > {path}{study}seg_M_phylo/NSmprot_align.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/Gnnucl.fasta > {path}{study}seg_M_phylo/Gnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/Gcnucl.fasta > {path}{study}seg_M_phylo/Gcnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}seg_M_phylo/NSmnucl.fasta > {path}{study}seg_M_phylo/NSmnucl_align.aln-fasta",shell=True)
