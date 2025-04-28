#Creation by CHARLES Seanna seanna.charles@univ-lyon1.fr
#Import
import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import AlignIO
from Bio import Phylo

import re

import subprocess

import random

import requests

from geopy.geocoders import Nominatim

import time

#The path may need to change
path = "/home/seannacharles/Documents/"

#######Function for the Script Genbank extract##################
#This part create a set of 100 random colors to create the itol legends 
cmap=plt.get_cmap("hsv")
Rcolor=[cmap(i/100)for i in range(100)] #Colors in rgb format
colors=['#%02x%02x%02x'%(int(r*255),int(g*255),int(b*255))for r,g,b,_ in Rcolor] #Colors in hexa format
random.shuffle(colors)


#This is the diffent type of legends
type_files=["Geolocalisation","Host","Serogroup"]

### Function Modif_ID 
# It is necessary to shorten the sequence identifier for dnadist (maximum 10 characters)
# It is important to uniformize the identifiers 
# It is important to remove all special characters that could cause problems in the future use of the software.
# It is important to remove terms that do not allow the different organisms to be differentiated (Add the terms you want in this function)

def Modif_ID(ID) : 
    #These changes are necessary to ensure the uniformity of the identifiers used by the various organisations.
    #Add you own modification if you feel the need
    ID=ID.lower()
    ID=ID.replace(" ","")
    ID=ID.replace("_","")
    ID=ID.replace("-","")
    ID=ID.replace("'","")
    ID=ID.replace("orthobunya","")
    ID=ID.replace("orthonairo","")
    ID=ID.replace("nairovirid","")
    ID=ID[:9]
    return ID

### Function Colour_choice
# A colour is assigned to each parameter that may appear in a legend
# The colours are chosen at random 
# All families will have the same captions
#Legencolors save for each legend label the associated colors to have the same color for the same label 
legencolors={}

def colour_choice (region):
    global colors
    if region not in legencolors : 
        choice=colors[0]
        colors.remove(choice)
        legencolors[region]=[choice]
        return(choice)
    else : 
        return(legencolors[region][0])

### Function modif_localisation
# For each geographical region, it is important to keep only those elements that identify the host's country of origin.

def modif_localisation(localisation) :
    localisation=localisation.split(":")[0]
    return localisation
### Function geographic_region
# For the captions, each country is replaced by its geographical region.
# At the moment I can't get pycountry_convert to work, which would simplify this function.
# If the geographical region can't be find the country will be kept

def geographic_region(pays) : 

    geolocator = Nominatim(user_agent=f"{pays}continentcharles.seanna@outlook.fr")
    location = geolocator.geocode(pays,timeout=10)


    latitude=location.latitude
    longitude=location.longitude

    if -35 <= latitude <= 37 and -20 <= longitude <= 55:
        return( "Africa")
    elif 35 <= latitude <= 70 and -10 <= longitude <= 40:
        return( "Europe")
    elif 5 <= latitude <= 55 and 70 <= longitude <= 140:
        return( "Asia")
    elif -55 <= latitude <= 15 and -80 <= longitude <= -35:
        return( "South America")
    elif 15 <= latitude <= 72 and -170 <= longitude <= -50:
        return( "North America")
    elif -50 <= latitude <= -10 and 110 <= longitude <= 180:
        return( "Oceania")
    else:
        return(pays)
   

### Function modif host
# For each host, it is important to keep only those elements that allow the species to be identified

def modif_host(host) : 
    if "culex" in host.lower() : 
        host ="Culex"
    elif "homo sapiens" in host.lower():
        host="Homo sapiens"
    elif "midge" in host.lower() :
        host="Midges"
    elif "tick" in host.lower() : 
        host="Tick"

    host.split(";")[0]
    host.replace("sp.",'')
    host.replace("sentinel.",'')
    host=re.sub(r"\([^)]*\)", "", host) #replace all the text betwenn ()

    return host

### Function deforder 
# For the captions, each species is replaced by its order
# If the order can't be find the species will be kept
# This function send a request to the NCBI database

def deforder(species) : 
    Entrez.email="charles.seanna@outlook.fr"

    try : 
        search= Entrez.esearch(db="taxonomy",term=species)
        record=Entrez.read(search)

        if not record["IdList"]:
            return species

        tax_id = record["IdList"][0]
        tax_info= Entrez.efetch(db="taxonomy",id=tax_id,retmode="xml")
        data=Entrez.read(tax_info)[0]

        for rank, name in zip(data["LineageEx"],data["LineageEx"]):
            if rank["Rank"]=="order" :
                return name["ScientificName"]

        return data["ScientificName"]

    except Exception as e : 
    
        return species

### Function write_file
# This function writes the caption files required for itol
# This function writes individual files for each nucleotide and protein sequence 

def write_file(Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,sp,serogroup,color_sero) : 
    global type_files
    global path

    organism.replace(" ","")
    name_geol=f"{path}{study}/{type_files[0]}_{study}_itol.txt"
    name_host=f"{path}{study}/{type_files[1]}_{study}_itol.txt"
    name_label=f"{path}{study}/Label_{study}_itol.txt"
    name_sero=f"{path}{study}/{type_files[2]}_{study}_itol.txt"

    geoloc=open(name_geol,'a')
    geoloc.write(f"##{localisation}\n")
    text=str (f"{ID}\trange\t{region_color}\t{region}\n")
    geoloc.write(text)
    geoloc.close()

    Host_f=open(name_host,'a')
    Host_f.write(f"##{host}\n")
    text=str (f"{ID}\trange\t{order_color}\t{order}\n")
    Host_f.write(text)
    Host_f.close()


    label=open(name_label,'a')
    label.write(f"##{organism}\n")
    text=str (f"{ID}\t{organism}\n")
    label.write(text)
    label.close()

    sero=open(name_sero,'a')
    sero.write(f"##{serogroup}\n")
    text=str (f"{ID}\trange\t{color_sero}\t{serogroup}\n")
    sero.write(text)
    sero.close()


    nom_prot_file=f"{strain}_{ID_prot}_{ID}"
    nom_prot_file.replace(" ","")
    nom_nucl_file=f"{strain}_{Nu}_{ID}"
    nom_nucl_file.replace(" ","")


    name_geol=f"{path}{study}/{type_files[0]}_{study}NSmGn_itol.txt"
    name_host=f"{path}{study}/{type_files[1]}_{study}NSmGn_itol.txt"
    name_label=f"{path}{study}/Label_{study}NSmGn_itol.txt"
    name_sero=f"{path}{study}/{type_files[2]}_{study}NSmGn_itol.txt"

    geoloc=open(name_geol,'a')
    geoloc.write(f"##{localisation}\n")
    text=str (f"N{ID}\trange\t{region_color}\t{region}\nG{ID}\trange\t{region_color}\t{region}\n")
    geoloc.write(text)
    geoloc.close()

    Host_f=open(name_host,'a')
    Host_f.write(f"##{host}\n")
    text=str (f"N{ID}\trange\t{order_color}\t{order}\nG{ID}\trange\t{order_color}\t{order}\n")
    Host_f.write(text)
    Host_f.close()


    label=open(name_label,'a')
    label.write(f"##{organism}\n")
    text=str (f"N{ID}\tN{organism}\nG{ID}\tG{organism}\n{ID}\t{organism}\n")
    label.write(text)
    label.close() 

    sero=open(name_sero,'a')
    sero.write(f"##{serogroup}\n")
    text=str (f"N{ID}\trange\t{color_sero}\t{serogroup}\nG{ID}\trange\t{color_sero}\t{serogroup}\n")
    sero.write(text)
    sero.close()

    nom_prot_file=f"{strain}_{ID_prot}_{ID}"
    nom_prot_file.replace(" ","")
    nom_nucl_file=f"{strain}_{Nu}_{ID}"
    nom_nucl_file.replace(" ","")







    if sp =="" : 
        with open (f"{path}{study}/M_nucl/{nom_nucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/M_prot/{nom_prot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")

    elif sp =="CCHFV" :
        with open (f"{path}{study}/M_nucl/{strain}_{Nu}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/M_prot/{strain}_{ID_prot}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")

    elif sp=="Nairoviridae" :
        with open (f"{path}{study}/M_nucl/{nom_nucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/M_prot/{nom_prot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")
            
    else :
        with open (f"{path}{study}/{nom_nucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/{nom_prot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")    

### Function first_call 
# This function writes the header for caption files for Itol
      
def first_call(study) : 
    global type_files
    global path
    for name in type_files : 
        name_files = f"{path}{study}/{name}_{study}_itol.txt"
        file = open(name_files,"w")
        file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATA\n")       
        file.close()
    for name in type_files : 
        name_files = f"{path}{study}/{name}_{study}NSmGn_itol.txt"
        file = open(name_files,"w")
        file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATA\n")       
        file.close()
    name_files = f"{path}{study}/Label_{study}_itol.txt"
    file = open(name_files,"w")
    file.write(f"LABELS\nSEPARATOR TAB\nDATA\n")       
    file.close()
    name_files = f"{path}{study}/Label_{study}NSmGn_itol.txt"
    file = open(name_files,"w")
    file.write(f"LABELS\nSEPARATOR TAB\nDATA\n")       
    file.close()

### Function linux_line 
# This function concatenates all the species for each study into one file 

def linux_line(sp,study) :
    subprocess.run(f"rm {path}{study}/M_prot_all_{sp}.fasta",shell=True)
    subprocess.run(f"rm {path}{study}/M_nucl_all_{sp}.fasta",shell=True)

    subprocess.run(f"cat {path}{study}/M_prot/* >> {path}{study}/M_prot_all_{sp}.fasta",shell=True)
    subprocess.run(f"cat {path}{study}/M_nucl/* >> {path}{study}/M_nucl_all_{sp}.fasta",shell=True)

### Function extract_id
# This function extracts the following information: Host, Host country, Strain, Nucleic identifier, Protein identifier for each genbank file.

name=f"{path}Metadata_studies.txt"
allmetadata=open(name,'w')
head=f"Virus name\tID itol tree\tTaxonomy\tSerogroup\tProtein accesion number\tSize Proteique sequence\tStart MLD\tEnd MLD\tSart GP38\tEnd GP38\tStart Gn\tEnd Gn\tStart NSm\tEnd NSm\tStart Gc\tEnd Gc\tNucleique accession number\tSize Nucelique sequence\tStrand\tStrain\tGeographical region\tHost\tHost source\tReferences\n"
allmetadata.write(head)
allmetadata.close()

def extract_id(Datas,study,sp) :

    global allmetadata

    for record in SeqIO.parse(Datas,"genbank"): 
        #Variable initialisation 
        host="Na"
        localisation="Na"
        strain="Na" 
        order="Na"
        region="Na"
        ID=''
        serogroup='Na'

        Nu= record.id
        Sequence_nu = record.seq
        sizenu=len(Sequence_nu)
        organism= record.annotations.get('organism',"Inconnu")
        organism=organism.replace('virus','')
        fullorganisme=organism

        taxonomy=record.annotations.get('taxonomy')[7:]
        taxonomy=' : '.join(taxonomy)

        for feature in record.features : 
            if "host" in feature.qualifiers : 
                host = feature.qualifiers['host'][0]
                host=modif_host(host)
                order = deforder(host) 

            if "country" in feature.qualifiers :
                localisation = feature.qualifiers['country'][0]
                localisation=modif_localisation(localisation)
                region = geographic_region(localisation)

            elif "geo_loc_name" in feature.qualifiers :
                localisation = feature.qualifiers['geo_loc_name'][0]
                localisation=modif_localisation(localisation)
                region = geographic_region(localisation)

            if "strain" in feature.qualifiers :
                strain = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['strain'][0])

            elif "isolate" in feature.qualifiers :
                strain = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['isolate'][0])

            if feature.type == "CDS" and "translation" in feature.qualifiers :
                Sequence_aa=feature.qualifiers["translation"][0] 
                sizeprot=len(Sequence_aa)    
                            
            if "protein_id" in feature.qualifiers :
                ID_prot = feature.qualifiers['protein_id'][0]

        output_prot=extract_prot_data(ID_prot,study)

        order_color=colour_choice(order)
        region_color=colour_choice(region)
        serogroup=find_serogroup(fullorganisme)
        color_sero=colour_choice(serogroup)

        if study=="CCHFV" : 
            ID=Modif_ID(str(ID_prot+strain))
            Description_nu = f"{ID}"
            Description_prot = f"{ID}"
            organism=ID

        else : 
            ID=Modif_ID(str(organism+ID_prot))
            Description_nu = f"{ID}"
            Description_prot = f"{ID}"
        
        Description_nu.replace(" ","")
        Description_prot.replace(" ","")



        write_file(Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,sp,serogroup,color_sero)
        
        
        line=f"{fullorganisme}\t{ID}\t{taxonomy}\t{serogroup}\t{ID_prot}\t{sizeprot}\t{output_prot}{Nu}\t{sizenu}\t\t{strain}\t{region}: {localisation}\t{order}: {host}\t\t\t\n"
        name=f"{path}Metadata_studies.txt"

        allmetadata=open(name,'a')
        allmetadata.write(line)
        allmetadata.close()
    return strain,ID_prot,Nu


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

def extract_prot_data(ID_prot,study) : 

    all_proteins=["MLD","GP38","Gn","NSm","Gc"]
    Entrez.email = "charles.seanna@outlook.fr"  
    handle = Entrez.efetch(db="protein", id=ID_prot, rettype="gp", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    time.sleep(2)  
    with open(f"{path}{study}/M_prot/{ID_prot}.gp", "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")
        output_handle

    record = SeqIO.read(f"{path}{study}/M_prot/{ID_prot}.gp", "genbank")
    proteins={}
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
        elif (prot=="GP38" or prot=="Gn") and "GP38+Gn" in proteins:
            output+=f'GP38+Gn: {str(proteins["GP38+Gn"][0])}\tGP38+Gn: {str(proteins["GP38+Gn"][1])}\t'
        else : 
            output+=f"\t\t"
    subprocess.run(f"rm {path}{study}/M_prot/{ID_prot}.gp",shell=True)
    return output

df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0)
df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)

def find_serogroup(organism) : 
    global df
    end_line=df.shape[0]
    for line in range(0,end_line) : 
        if df.iloc[line,0].lower() in organism.lower() or organism.lower() in df.iloc[line,0].lower() : 
            return str(df.iloc[line,1])


#######Function for the Script mark_prot##################


##Function compt_ali_nucl
# This function uses the position of the proteins of the target organism (proteic sequence) to identify the position of the proteins in the alignment.
# It is important not to just translate the protein sequence into RNA in order to keep the intergenic regions.

def compt_ali_nucl(sequence,debseq,pos_prot) :
    pos_align=[]
    compt_ali=0
    compt_nu=0
    compt_prot=0
    compt_codon=0
    codon=""
    deb=True
    for nu in sequence: 
        compt_ali+=1
        if nu!= "-" : 
            compt_nu+=1
            if compt_nu >=debseq: 
                compt_codon+=1
                codon+=nu
                if compt_codon==3 : 
                    compt_codon=0
                    compt_prot+=1
                    if compt_prot in pos_prot :
                            
                        if deb==True : 
                            deb=False
                            pos_align.append(compt_ali-2)
                            pos_prot.remove(compt_prot)
                        else : 
                            deb=True
                            pos_align.append(compt_ali)
                            pos_prot.remove(compt_prot) 
                        
                    codon=""
    return pos_align
  
##Function compt_alignement
# This function uses the position of the target organism's proteins (sequence proteic) to identify the position of the proteins in the alignment.

def compt_alignement (sequence,pos_prot) : 
    
    pos_align=[]
    compt_ali=0
    compt_prot=0
    mark_prot=""
    deb=True
    for aa in sequence : 
        if aa=='-' : 
            compt_ali+=1
        else : 
            compt_ali+=1
            compt_prot+=1

        if compt_prot in pos_prot :
            if deb==True : 
                mark_prot+="|"
                mark_prot+=aa
                deb=False
                pos_align.append(compt_ali)
                pos_prot.remove(compt_prot)
            else : 
                mark_prot+=aa
                mark_prot+="|"
                pos_align.append(compt_ali)
                deb=True
                pos_prot.remove(compt_prot)
        else : 
            mark_prot+=aa
    return pos_align

##Function mark_sequence
#This function marks the start and end of proteins in protein and nucleic acid sequences. 

def mark_sequence(sequence,pos_align,name,description) : 
    compt_ali=0
    mark_prot=""
    deb=True
    for aa in sequence : 
        compt_ali+=1
        if compt_ali in pos_align :
            if deb==True : 
                mark_prot+="|"
                mark_prot+=aa
                deb=False
             
            else : 
                mark_prot+=aa
                mark_prot+="|"
                deb=True
        else : 
            mark_prot+=aa
    
    file = open(name,"a")
    text=str (">"+description +'\n'+mark_prot+'\n')
    file.write(text)
    file.close()




#######Function for the split_protein##################


##Function Extract_data
# This function separates the different proteins that have been tagged for all sequences.

def Extract_data(file,sp,Protname,study) : 
    if "nucl" in sp  :
        AA_L=["a","t","c","g","-","n"]
    else :
        AA_L=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","-","|"]
    Data={}

    for record in SeqIO.parse(f"{path}{study}/{file}","fasta") :
        if record.id not in Data :
            seq=re.findall(r'\|(.*?)\|',str (record.seq))
            for index,prot in enumerate(seq) : 
                for AA in prot:
                    if AA not in AA_L : 
                        prot_mod=prot.replace(str(AA),"n")
                        seq[index]=prot_mod
                    
            Data[record.description]=seq
    write_file_prot(Data,sp,Protname,study)


## Function write_file_prot
# This function writes a file for each protein with all the sequences studied.  
 
def write_file_prot(Data,sp,Protname,study) : 
 
    for ID,seq in Data.items() : 
        for index,prot in enumerate(seq) : 
            name=f"{path}{study}/seg_M_phylo/"+str(Protname[index])+sp+".fasta"
            file = open(name,"a")
            text=str (">"+ID+'\n'+Seq(prot)+'\n')
            file.write(text)
            file.close()
            
#######Function for the Script delete_dbl##################
##Function select_doublon
#This function identified all the possible duplicate sequences
#The possible duplicates are identified using a specific threshold; this threshold was found using the dnadist file and the visible rake of the phylogenetic tree (specific to each study and not provided by this script). 

def select_doublon(protein,study,target,ths) :
    global path
    df = pd.read_csv(f"{path}{study}/seg_M_phylo/{protein}nucl_dnadist.txt", delimiter="\t",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(float)
    
    
    Doublon=[]
    end_line=df.shape[0]
    end_col=df.shape[1]
    line=0
    col=0
    while line <end_line:
        col=line+1
        while col <end_col-1 :
            if df.iloc[line,col] <=ths and df.iloc[line,col]!=-1  :
                if target not in df.index[line]:
                    ID=str(df.index[line])
                    Doublon.append(ID)
                    df=df.drop(columns=ID)
                    df=df.drop(index=ID)
                    #df = df.reset_index(drop=True)
                    end_line=df.shape[0]
                    end_col=df.shape[1]
                    line=0

            col+=1
            

        line+=1
    noms_colonnes = df.columns
    return Doublon,noms_colonnes
    #delete_dbl(Doublon,study,protein) 

##Function col_name
#This function delete the duplicate sequence name 
    
def col_name(old_col_names,col_name,Protein) :

    if len(old_col_names) >0 :
        if Protein=="NSmGn" : 
            col_name=[element[-9:] for element in col_name]

            old_col_names = [element for element in old_col_names if element in col_name]
        else :
            old_col_names = [element for element in old_col_names if element in col_name]
    else : 
        old_col_names=col_name

    return old_col_names

##Function delete_dbl
#This function delete the identified duplicates for each protein for each study

def delete_dbl(Doublon,study,protein,old_col_names):
    global path
    
    
    name=f"{path}{study}/seg_M_phylo/{protein}nucl_sdbl.fasta"
    
    subprocess.run(f"rm {name}",shell=True)

    if study!="CCHFV" :
        file_nucl=f"{path}{study}/seg_M_phylo/{protein}nucl_trimal.fasta"
        file_prot=f"{path}{study}/seg_M_phylo/{protein}prot_trimal.fasta"
    else : 
        file_nucl=f"{path}{study}/seg_M_phylo/{protein}nucl.aln-fasta"
        file_prot=f"{path}{study}/seg_M_phylo/{protein}prot.aln-fasta"

    for record in SeqIO.parse(f"{file_nucl}","fasta") :
        
        if record.description not in Doublon and record.description in old_col_names :
            file = open(name,"a")
            text=str (">"+record.description +'\n'+record.seq+'\n')
            file.write(text)
            file.close()

    output_phylip=f"{path}{study}/seg_M_phylo/{protein}nucl_sdbl.phy"
    subprocess.run(f" mafft --auto {name} > {path}{study}/seg_M_phylo/{protein}nucl_sdbl.aln-fasta ",shell=True)
    alignment = AlignIO.read(f"{path}{study}/seg_M_phylo/{protein}nucl_sdbl.aln-fasta", "fasta")
    AlignIO.write(alignment, output_phylip, "phylip")


    name=f"{path}{study}/seg_M_phylo/{protein}prot_sdbl.fasta"
    subprocess.run(f" rm {name}",shell=True)
    
    for record in SeqIO.parse(f"{file_prot}","fasta") :
        
        if record.description not in Doublon and record.description in old_col_names :
            file = open(name,"a")
            text=str (">"+record.description+'\n'+record.seq+'\n')
            file.write(text)
            file.close()

    output_phylip=f"{path}{study}/seg_M_phylo/{protein}prot_sdbl.phy"
    subprocess.run(f" mafft --auto {name} > {path}{study}/seg_M_phylo/{protein}prot_sdbl.aln-fasta ",shell=True)
    alignment = AlignIO.read(f"{path}{study}/seg_M_phylo/{protein}prot_sdbl.aln-fasta", "fasta")
    AlignIO.write(alignment, output_phylip, "phylip")
    

