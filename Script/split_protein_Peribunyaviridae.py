import subprocess
import sys

from All_Functions import Extract_data, write_file

path =str(sys.argv[1])
study= str(sys.argv[2])
res_file= str(sys.argv[3])


 

print('split prot Peribunyaviridae')


Protname=["Gn","NSm","Gc"]
Extract_data("prot_all_align_marke_prot.aln-fasta","protSH",Protname,study)
Extract_data("nucl_all_align_marke.aln-fasta","nuclSH",Protname,study)

Extract_data("prot_all_align_marke_prot_hedwig.fasta","prot",Protname,study)
Extract_data("nucl_all_align_marke_hedwig.fasta","nucl",Protname,study)

Protname=["Gn","Gc"]
Extract_data("prot_marke_prot_khurdun.fasta","protKhurdun",Protname,study)
Extract_data("nucl_marke_khurdun.fasta","nuclKhurdun",Protname,study)






subprocess.run(f"mafft --auto {path}{study}{res_file}Gnprot.fasta > {path}{study}{res_file}Gnprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}Gcprot.fasta > {path}{study}{res_file}Gcprot_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}NSmprot.fasta > {path}{study}{res_file}NSmprot_align.aln-fasta",shell=True)


subprocess.run(f"mafft --auto {path}{study}{res_file}Gnnucl.fasta > {path}{study}{res_file}Gnnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}Gcnucl.fasta > {path}{study}{res_file}Gcnucl_align.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}{res_file}NSmnucl.fasta > {path}{study}{res_file}NSmnucl_align.aln-fasta",shell=True)
