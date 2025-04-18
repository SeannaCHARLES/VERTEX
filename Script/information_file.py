from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import AlignIO
from Bio import Phylo
import re
import pandas as pd

path="/home/seannacharles/Documents/"
studies=["Perybunyaviridae"]
name=f"{path}Metadata_studies.txt"
metadata=open(name,'w')
head=f"Virus name\tAbreviation\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteique sequence\tStart Gn\tEnd Gn\tStart Gc\tEnd Gc\tStart NSm\tEnd NSm\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tNucleique accession number\tSize Nucelique sequence\tStart Gn\tEnd Gn\tStart Gc\tEnd Gc\tStart NSm\tEnd NSm\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tStrand\tStrain\tGeographical region\tHost\tHost source\tReferences\n"
metadata.write(head)
metadata.close()


from Bio import Entrez, SeqIO

Entrez.email = "charles.seanna@outlook.fr"  



all_proteins=["MLD","GP38","Gn","NSm","Gc"]

def extract_prot_data(ID_prot) : 
    handle = Entrez.efetch(db="protein", id=ID_prot, rettype="gp", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    with open(f"{ID_prot}.gp", "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")
        output_handle

    record = SeqIO.read(f"{ID_prot}.gp", "genbank")
    proteins={}
    find_serogroup(record.annotations.get('organism',"Inconnu"))
    for feature in record.features :

        if 'region_name' in feature.qualifiers : 
            prot=(feature.qualifiers['region_name'][0])
            start=(feature.location.start)
            end=(feature.location.end)
            prot=modif_name_prot(prot) 
            proteins[prot]=[start,end]    
    output=""        
    for prot in all_proteins : 
        if prot in proteins : 
            output+=f'{str(proteins[prot][0])}\t{str(proteins[prot][1])}\t'
        elif prot=="GP38" or prot=="Gn" and "GP38+Gn" in proteins.keys():
            output+=f'GP38+Gn: {str(proteins["GP38+Gn"][0])}\tGP38+Gn: {str(proteins["GP38+Gn"][1])}\t'
        else : 
            output+=f"\t\t"
    return output

def modif_name_prot(name_prot) : 
    name_prot=name_prot.upper()
    if "HANTA_G2" in name_prot  : 
        name_prot="Gc"
    name_prot=name_prot.replace("BUNYA","")
    name_prot=name_prot.replace("_","")
    name_prot=name_prot.replace("G2","Gn")
    name_prot=name_prot.replace("G1","Gc")
    name_prot=name_prot.replace("NSM","NSm")
    name_prot=name_prot.replace("DUF5585","MLD")
    name_prot=name_prot.replace("NAIROVIRUSM","GP38+Gn")
    return name_prot

df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0)
df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)

def find_serogroup(organism) : 
    global df
    end_line=df.shape[0]
    for line in range(1,end_line) : 
        if df.iloc[line,1].lower() in organism.lower() : 
            print(f'{df.iloc[line,0]}\t{df.iloc[line,1]}')
            return str(df.iloc[line,1])
output=extract_prot_data("XCY99298.1")





