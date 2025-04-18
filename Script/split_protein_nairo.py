from Bio import SeqIO
import re
import subprocess
from Bio.Seq import Seq
from Bio import AlignIO
from All_Function import Extract_data, write_file
path = "/home/seannacharles/Documents/Nairoviridae/"
study="Nairoviridae"

Protname=["MLD","GP38","Gn","NSm","Gc"]


Extract_data("M_prot_all_align_marke_prot.fasta","prot",Protname,study)
Extract_data("M_nucl_all_align_marke.fasta","nucl",Protname,study)


subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnprot.fasta > {path}seg_M_phylo/Gnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcprot.fasta > {path}seg_M_phylo/Gcprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmprot.fasta > {path}seg_M_phylo/NSmprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38prot.fasta > {path}seg_M_phylo/GP38prot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDprot.fasta > {path}seg_M_phylo/MLDprot_align.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnnucl.fasta > {path}seg_M_phylo/Gnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcnucl.fasta > {path}seg_M_phylo/Gcnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmnucl.fasta > {path}seg_M_phylo/NSmnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38nucl.fasta > {path}seg_M_phylo/GP38nucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDnucl.fasta > {path}seg_M_phylo/MLDnucl_align.aln-fasta",shell=True)

