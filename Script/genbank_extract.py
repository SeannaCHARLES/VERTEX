
import subprocess
import sys 

from All_Functions import first_call,concatenates,extract_id

path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
studiesspeOG=list(sys.argv[3].split(" "))
names=list(sys.argv[4].split(" "))
NUIDspeOG=list(sys.argv[5].split(" "))
email=sys.argv[6]
list_prot = [x.split(';') for x in sys.argv[7].split()]
list_prot_speOG = [x.split(';') for x in sys.argv[8].split()]
positionspeOG=list(sys.argv[8].split(" "))
OG_prot=list(sys.argv[9].split(" "))
Diff_name=list(sys.argv[10].split(";"))
chain_prot=""

# print(email)
# print(positionspeOG)
# print(list_prot_speOG)
# print(Diff_name)
vv_prot={}
for study,prot in zip(studies,list_prot) :
    print(study)
    subprocess.run(f"rm {path}{study}/nucl/* {path}{study}/prot/* >/dev/null",shell=True)
    chain_prot=""
    vv_prot[study]=prot

    for p in prot : 

        chain_prot+=f"Start{p}\tEnd{p}\t"
    subprocess.run(f"echo 'Study\tVirus name\tID itol tree\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteic sequence\t{chain_prot}Nucleic accession number\tSize Nucleic sequence\tStrand\tStrain\tGeographical region\tHost or vector\tReferences' >{path}{study}/Metadata_{study}.txt",shell=True)

    subprocess.run(f"echo 'Study\tID itol tree\t{chain_prot}' >{path}{study}/Pos_prot_{study}.txt",shell=True)
    with open(f"{path}{study}/sequences_{study}.gb","r") as Datas:
        sp=''
        first_call(path,study)
        extract_id(path,Datas,study,sp,prot,email,Diff_name,chain_prot[:-1].split("\t"),vv_prot[study])
        concatenates(path,study)

    subprocess.run(f"mafft --auto {path}{study}/prot_all_{study}.fasta > {path}{study}/prot_all_align_{study}.aln-fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/nucl_all_{study}.fasta > {path}{study}/nucl_all_align_{study}.aln-fasta",shell=True)

for study,ID,prot in zip(studiesspeOG, NUIDspeOG,list_prot_speOG):
    
    with open(f"{path}{study}/{ID}.gb","r") as Datas:
        sp='OGspe'
        extract_id(path,Datas,study,sp,prot,email,Diff_name,chain_prot[:-1].split("\t"),vv_prot[study])
