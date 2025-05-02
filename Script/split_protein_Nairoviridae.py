import subprocess
import sys 

from All_Functions import Extract_data, write_file

path =str(sys.argv[1])
study= str(sys.argv[2])
res_file= str(sys.argv[3])

Protname=["MLD","GP38","Gn","NSm","Gc"]


Extract_data("prot_all_align_marke_prot.fasta","prot",Protname,study)
Extract_data("nucl_all_align_marke.fasta","nucl",Protname,study)


subprocess.run(f"mafft --auto {path}{study}{res_file}Gnprot.fasta > {path}{study}{res_file}Gnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}Gcprot.fasta > {path}{study}{res_file}Gcprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}NSmprot.fasta > {path}{study}{res_file}NSmprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}GP38prot.fasta > {path}{study}{res_file}GP38prot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}MLDprot.fasta > {path}{study}{res_file}MLDprot_align.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}{study}{res_file}Gnnucl.fasta > {path}{study}{res_file}Gnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}Gcnucl.fasta > {path}{study}{res_file}Gcnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}NSmnucl.fasta > {path}{study}{res_file}NSmnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}GP38nucl.fasta > {path}{study}{res_file}GP38nucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}MLDnucl.fasta > {path}{study}{res_file}MLDnucl_align.aln-fasta",shell=True)

