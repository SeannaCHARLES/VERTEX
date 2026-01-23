
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
Diff_name=list(sys.argv[10].split(" "))
Diff_name_OG=list(sys.argv[11].split(" "))
speci_order=list(sys.argv[12].split(";"))
orgaWIWO=list(sys.argv[13].replace(" ", "\n").splitlines())
INOUT=list(sys.argv[14].replace(" ", "\n").splitlines())
protINOUT=list(sys.argv[15].replace(" ", "\n").splitlines())
Allstudies=list(sys.argv[16].split(" "))
ALLNUID=list(sys.argv[17].split(" "))
namespeprot=list(sys.argv[18].split(" "))
namespeprotOG=list(sys.argv[19].split(" "))
for pos,element in enumerate(namespeprot) : 
    namespeprot[pos]=element.split(";")
for pos,element in enumerate(namespeprotOG) : 
    namespeprotOG[pos]=element.split(";")

chain_prot=""

dicy={}
for all_s, IO,prots,virus,list_ID in zip(Allstudies,INOUT,protINOUT,orgaWIWO,ALLNUID) :
   
    list_ID=list_ID.split(";")
    #metadata = pd.read_csv(f"{path}{all_s}/Metadata_{all_s}.txt", delimiter="\t",header=0,index_col=0)
    
    if IO!='' :
        if all_s not in dicy : 
            dicy[all_s]={}
       
        if IO=="OUT" : 
            for v in virus.split(";") : 
                if v not in dicy[all_s] :
                    dicy[all_s][v]=[prots] 
                else : 
                    dicy[all_s][v].append([prots])

        elif IO=="IN" : 
            for  ID in list_ID :
                if ID not in virus.split(";") :
                    if ID not in dicy[all_s] :
                        dicy[all_s][ID]=[prots]
                    elif prots not in  dicy[all_s][ID] : 
                        dicy[all_s][ID].append(prots)
        elif IO=='none' :
            if all_s not in dicy : 
                dicy[all_s]={}

if not dicy : 
    for all_s in Allstudies : 
        dicy[all_s]=["-","-","-"]



vv_prot={}

for study,prot,prot_names,spe_prot in zip(studies,list_prot,Diff_name,namespeprot) :
    
    subprocess.run(f"rm {path}{study}/nucl/* {path}{study}/prot/* >/dev/null",shell=True)
    chain_prot=""
    vv_prot[study]=prot

    for p in prot : 

        chain_prot+=f"Start{p}\tEnd{p}\t"
    subprocess.run(f"echo 'Study\tVirus name\tID itol tree\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteic sequence\t{chain_prot}Nucleic accession number\tSize Nucleic sequence\tStrand\tStrain\tGeographical region\tHost or vector\tReferences' >{path}{study}/Metadata_{study}.txt",shell=True)

    subprocess.run(f"echo 'Study\tID itol tree\t{chain_prot}' >{path}{study}/Pos_prot_{study}.txt",shell=True)
    with open(f"{path}{study}/sequences_{study}.gb","r") as Datas:
        sp=''
        first_call(path,study,spe_prot)
        
        extract_id(path,Datas,study,sp,prot,email,prot_names.split(";"),vv_prot[study],speci_order,dicy,spe_prot)
        concatenates(path,study)

    subprocess.run(f"mafft --auto {path}{study}/prot_all_{study}.fasta > {path}{study}/prot_all_align_{study}.aln-fasta",shell=True)
    subprocess.run(f"mafft --auto {path}{study}/nucl_all_{study}.fasta > {path}{study}/nucl_all_align_{study}.aln-fasta",shell=True)

for study,ID,prot,prot_names,spe_prot in zip(studiesspeOG, NUIDspeOG,list_prot_speOG,Diff_name_OG,namespeprotOG):
    
    with open(f"{path}{study}/{ID}.gb","r") as Datas:
        sp='OGspe'
        extract_id(path,Datas,study,sp,prot,email,prot_names.split(";"),vv_prot[study],speci_order,dicy,spe_prot)
