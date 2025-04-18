# Script qui sépare les différentes proteines du segement M

import subprocess

from All_Function import Extract_data, write_file

print('split prot CCHFV')
path = "/home/seannacharles/Documents/CCHFV/"
study="CCHFV"

Protname=["MLD","GP38","Gn","NSm","Gc"]


Extract_data(f"M_prot_all_align_CCHFV_marke_prot.aln-fasta","prot",Protname,study)

Extract_data(f"M_nucl_all_align_CCHFV_marke.aln-fasta","nucl",Protname,study)



subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnprot.fasta > {path}seg_M_phylo/Gnprot.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcprot.fasta > {path}seg_M_phylo/Gcprot.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmprot.fasta > {path}seg_M_phylo/NSmprot.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38prot.fasta > {path}seg_M_phylo/GP38prot.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDprot.fasta > {path}seg_M_phylo/MLDprot.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}seg_M_phylo/Gnnucl.fasta > {path}seg_M_phylo/Gnnucl.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/Gcnucl.fasta > {path}seg_M_phylo/Gcnucl.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/NSmnucl.fasta > {path}seg_M_phylo/NSmnucl.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/GP38nucl.fasta > {path}seg_M_phylo/GP38nucl.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}seg_M_phylo/MLDnucl.fasta > {path}seg_M_phylo/MLDnucl.aln-fasta",shell=True)

