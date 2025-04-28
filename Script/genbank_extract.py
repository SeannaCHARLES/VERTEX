
import subprocess

from All_Functions import Modif_ID,colour_choice,modif_localisation,geographic_region,modif_host,deforder,write_file,first_call,linux_line,extract_id

path = "/home/seannacharles/Documents/"

study= "CCHFV"

with open(f"{path}{study}/sequences_{study}_segM.gb","r") as Datas:
    
    sp=f"{study}"

    first_call(study)
    extract_id(Datas,study,sp)
    linux_line(sp,study)

subprocess.run(f"mafft --auto {path}{study}/M_prot_all_{sp}.fasta > {path}{study}/M_prot_all_align_{sp}.aln-fasta",shell=True)
subprocess.run(f"mafft --auto {path}{study}/M_nucl_all_{sp}.fasta > {path}{study}/M_nucl_all_align_{sp}.aln-fasta",shell=True)


study= "Perybunyaviridae"


with open(f"{path}{study}/sequences_{study}_segM.gb","r") as Datas:
    
    sp=""

    first_call(study)
    extract_id(Datas,study,sp)
    linux_line(sp,study)

    subprocess.run(f"mafft --auto {path}{study}/M_prot_all_.fasta > {path}{study}/M_prot_all_align_{sp}.aln-fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/M_nucl_all_.fasta > {path}{study}/M_nucl_all_align_{sp}.aln-fasta",shell=True)


with open(f"{path}{study}/OV054446.1.gb","r") as Datas:
    
    sp="hedwig"

    extract_id(Datas,study,sp)

    subprocess.run(f"cat {path}{study}/M_prot_all_.fasta {path}{study}/ED-I-079-18_CAH0473050.1_hedwigcah.fasta >{path}{study}/M_prot_all_align_{sp}.fasta",shell=True)
    subprocess.run(f"cat {path}{study}/M_nucl_all_.fasta {path}{study}/ED-I-079-18_OV054446.1_hedwigcah.fasta >{path}{study}/M_nucl_all_align_{sp}.fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/M_prot_all_align_{sp}.fasta > {path}{study}/M_prot_all_align_{sp}.aln-fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/M_nucl_all_align_{sp}.fasta > {path}{study}/M_nucl_all_align_{sp}.aln-fasta",shell=True)

with open(f"{path}{study}/NC_078560.1.gb","r") as Datas:
    
    sp="kurdun"
    extract_id(Datas,study,sp)





study= "Nairoviridae"
print(f"Data extraction from genbank files for {study}")
with open(f"{path}{study}/sequences_{study}_segM.gb","r") as Datas:
    sp=f"{study}"

    first_call(study)
    extract_id(Datas,study,sp)
    linux_line(sp,study)

    subprocess.run(f"mafft --auto {path}{study}/M_prot_all_{sp}.fasta > {path}{study}/M_prot_all_align_{sp}.aln-fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/M_nucl_all_{sp}.fasta > {path}{study}/M_nucl_all_align_{sp}.aln-fasta",shell=True)


