
import subprocess
import sys 

from All_Functions import Modif_ID,colour_choice,modif_localisation,geographic_region,modif_host,deforder,write_file,first_call,linux_line,extract_id

path = sys.argv[1]
#path = "/home/seannacharles/"
studies=["CCHFV","Nairoviridae","Peribunyaviridae","hedwig","khurdun"]
subprocess.run(f"echo 'Virus name\tID itol tree\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteic sequence\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tStart Gn\tEnd Gn\tStart NSm\tEnd NSm\tStart Gc\tEnd Gc\tNucleic accession number\tSize Nucleic sequence\tStrand\tStrain\tGeographical region\tHost or vector\tReferences' >{path}Metadata_studies.txt",shell=True)
subprocess.run(f"echo 'ID itol tree\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tStart Gn\tEnd Gn\tStart NSm\tEnd NSm\tStart Gc\tEnd Gc\t' >{path}Prot_file_size.txt",shell=True)

for study in studies :
    sp= study
    print(f"Data extraction from genbank files for {study}")
    if study=="Peribunyaviridae" :
        with open(f"{path}{study}/sequences_{study}.gb","r") as Datas:   
            sp=""

            first_call(study)
            extract_id(Datas,study,sp)
            linux_line(sp,study)

            subprocess.run(f"mafft --auto {path}{study}/prot_all_.fasta > {path}{study}/prot_all_align_{sp}.aln-fasta",shell=True)
            subprocess.run(f"mafft --auto {path}{study}/nucl_all_.fasta > {path}{study}/nucl_all_align_{sp}.aln-fasta",shell=True)

    elif study=="hedwig" :
        study="Peribunyaviridae"
        with open(f"{path}{study}/OV054446.1.gb","r") as Datas:

            extract_id(Datas,study,sp)

            subprocess.run(f"cat {path}{study}/prot_all_.fasta {path}{study}/ED-I-079-18_CAH0473050.1_hedwigcah.fasta >{path}{study}/prot_all_align_{sp}.fasta",shell=True)
            subprocess.run(f"cat {path}{study}/nucl_all_.fasta {path}{study}/ED-I-079-18_OV054446.1_hedwigcah.fasta >{path}{study}/nucl_all_align_{sp}.fasta",shell=True)
            subprocess.run(f"mafft --auto {path}{study}/prot_all_align_{sp}.fasta > {path}{study}/prot_all_align_{sp}.aln-fasta",shell=True)
            subprocess.run(f"mafft --auto {path}{study}/nucl_all_align_{sp}.fasta > {path}{study}/nucl_all_align_{sp}.aln-fasta",shell=True)

    elif study=="khurdun" : 
        study="Peribunyaviridae"
        with open(f"{path}{study}/NC_078560.1.gb","r") as Datas:

            extract_id(Datas,study,sp)
    else : 
        with open(f"{path}{study}/sequences_{study}.gb","r") as Datas:

            sp=f"{study}"

            first_call(study)
            extract_id(Datas,study,sp)
            linux_line(sp,study)

        subprocess.run(f"mafft --auto {path}{study}/prot_all_{sp}.fasta > {path}{study}/prot_all_align_{sp}.aln-fasta",shell=True)
        subprocess.run(f"mafft --auto {path}{study}/nucl_all_{sp}.fasta > {path}{study}/nucl_all_align_{sp}.aln-fasta",shell=True)
    