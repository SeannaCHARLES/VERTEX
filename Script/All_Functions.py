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
#You email addresse for NCBI request
mail="charles.seanna@outlook.fr"

#The path may need to change
path = "/home/seannacharles/"
res_file="/res_phylo/"
#######Function for the Script Genbank extract##################
#This part create a set of 100 random colors to create the itol legends 
cmap=plt.get_cmap("hsv")
Rcolor=[cmap(i/60)for i in range(60)] #Colors in rgb format
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
legencolors["Na"]=["#FFFFFF"]
name_legend=f"{path}/All_legend_caption.txt"
All_legend=open(name_legend,'w')
text=str (f"Legend\tColor\nNa\t#FFFFFF\n")
All_legend.write(text)
All_legend.close()

def colour_choice (caption):
    global colors
    global path
    if caption not in legencolors : 
        choice=colors[0]
        colors.remove(choice)
        legencolors[caption]=[choice]

        name_legend=f"{path}/All_legend_caption.txt"
        All_legend=open(name_legend,'a')
        text=str (f"{caption}\t{choice}\n")
        All_legend.write(text)
        All_legend.close()
        return(choice)
    else : 
        return(legencolors[caption][0])

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
        return "Africa"
    elif 35 <= latitude <= 70 and -10 <= longitude <= 40:
        return "Europe"
    elif 35 <= latitude <= 75 and 40 < longitude <= 70:
        return "Eurasia"  
    elif 50 <= latitude <= 75 and 40 <= longitude <= 180:
        return "Russia"
    elif 5 <= latitude <= 55 and 70 <= longitude <= 140:
        return "Asia"
    elif -55 <= latitude <= 15 and -80 <= longitude <= -35:
        return "South America"
    elif 15 <= latitude <= 72 and -170 <= longitude <= -50:
        return "North America"
    elif -50 <= latitude <= -10 and 110 <= longitude <= 180:
        return "Oceania"
    else:
        return pays
   
### Function modif host
# For each host, it is important to keep only those elements that allow the species to be identified

def modif_host(host) : 
    if "culex" in host.lower() : 
        host ="Culex"
    elif "homo sapiens" in host.lower():
        host="Homo sapiens"
    elif "midge" in host.lower() :
        host="Midge"
    elif "tick" in host.lower() : 
        host="Ticks"
    elif "pigeon" in host.lower():
        host="Columba livia"
    elif "hamster" in host.lower():
        host="Cricetinae"
    
    host=host.split(";")[0]
    host=host.replace("sp.",'')
    host=host.replace("sentinel.",'')
    host=re.sub(r"\([^)]*\)", "", host) #replace all the text betwenn ()

    return host

### Function deforder 
# For the captions, each species is replaced by its order
# If the order can't be find the species will be kept
# This function send a request to the NCBI database

def deforder(species) : 
    global mail
    Entrez.email=mail
    sp =species
    try : 
        search= Entrez.esearch(db="taxonomy",term=species)
        record=Entrez.read(search)
        
        if not record["IdList"]:
            if sp not in species:
                return sp
            elif (sp!="" or sp!=" " or sp==species)and len(species.split()) > 1 :
                for part in species.split(" ") : 
                    spepart=deforder(part)
                    if spepart not in species : 
                        return spepart
                print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add a If statement in the modif_host function such as.\n elif {species.lower()} in host.lower() : \n\thost='Scientificname'")
                return species
            else  : 
                print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add a If statement in the modif_host function such as.\n elif {species.lower()} in host.lower() : \n\thost='Scientificname'")
                return species

        tax_id = record["IdList"][0]
        tax_info= Entrez.efetch(db="taxonomy",id=tax_id,retmode="xml")
        data=Entrez.read(tax_info)[0]

        for rank, name in zip(data["LineageEx"],data["LineageEx"]):
            if rank["Rank"]=="order" :
                return name["ScientificName"]

        return data["ScientificName"]

    except Exception as e : 
        if sp not in species:
            return sp
        elif (sp!="" or sp!=" " or sp==species)and len(species.split()) > 1 :
            for part in species.split(" ") : 
                spepart=deforder(part)
                if spepart not in species : 
                    return spepart
            print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add a If statement in the modif_host function such as.\n If {species.lower()} in host.lower() : \n\thost='Scientificname'")
            return species
        else  : 
            print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add a If statement in the modif_host function such as.\n If {species.lower()} in host.lower() : \n\thost='Scientificname'")
            return species

### Function write_file
# This function writes the caption files required for itol
# This function writes individual files for each nucleotide and protein sequence 

def write_file(Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,sp,serogroup,color_sero) : 
    global type_files
    global path

    organism.replace(" ","")
    name_geol=f"{path}{study}/legend/{type_files[0]}_{study}_itol.txt"
    name_host=f"{path}{study}/legend/{type_files[1]}_{study}_itol.txt"
    name_label=f"{path}{study}/legend/Label_{study}_itol.txt"
    name_sero=f"{path}{study}/legend/{type_files[2]}_{study}_itol.txt"
    name_host_images=f"{path}{study}/legend/Host_images{study}_itol.txt"

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

    HostI_f=open(name_host_images,'a')
    HostI_f.write(f"##{host}\n")
    text=str (f"{ID}\t-1\t1\t0\t0\t0\thttps://raw.githubusercontent.com/SeannaCHARLES/Images/main/{order}.png\n")
    HostI_f.write(text)
    HostI_f.close()

    nonucl_file=f"{strain}_{ID_prot}_{ID}"
    nonucl_file.replace(" ","")
    nonucl_file=f"{strain}_{Nu}_{ID}"
    nonucl_file.replace(" ","")


    name_geol=f"{path}{study}/legend/NSmGn/{type_files[0]}_{study}NSmGn_itol.txt"
    name_host=f"{path}{study}/legend/NSmGn/{type_files[1]}_{study}NSmGn_itol.txt"
    name_label=f"{path}{study}/legend/NSmGn/Label_{study}NSmGn_itol.txt"
    name_sero=f"{path}{study}/legend/NSmGn/{type_files[2]}_{study}NSmGn_itol.txt"
    name_host_images=f"{path}{study}/legend/NSmGn/Host_images{study}NSmGn_itol.txt"

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
    text=str (f"N{ID}\tN{organism}\nG{ID}\tG{organism}\n")
    label.write(text)
    label.close() 

    sero=open(name_sero,'a')
    sero.write(f"##{serogroup}\n")
    text=str (f"N{ID}\trange\t{color_sero}\t{serogroup}\nG{ID}\trange\t{color_sero}\t{serogroup}\n")
    sero.write(text)
    sero.close()

    HostI_f=open(name_host_images,'a')
    HostI_f.write(f"##{host}\n")
    text=str (f"N{ID}\t-1\t1\t0\t0\t0\thttps://raw.githubusercontent.com/SeannaCHARLES/Images/main/{order}.png\nG{ID}\t-1\t1\t0\t0\t0\thttps://raw.githubusercontent.com/SeannaCHARLES/Images/main/{order}.png\n")
    HostI_f.write(text)
    HostI_f.close()


    nomprot_file=f"{strain}_{ID_prot}_{ID}"
    nomprot_file=nomprot_file.replace(" ","")
    nomnucl_file=f"{strain}_{Nu}_{ID}"
    nomnucl_file=nomnucl_file.replace(" ","")







    if sp =="" : 
        with open (f"{path}{study}/nucl/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/prot/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")

    elif sp =="CCHFV" :
        with open (f"{path}{study}/nucl/{strain}_{Nu}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/prot/{strain}_{ID_prot}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")

    elif sp=="Nairoviridae" :
        with open (f"{path}{study}/nucl/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/prot/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")
            
    else :
        with open (f"{path}{study}/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")    

### Function first_call 
# This function writes the header for caption files for Itol
      
def first_call(study) : 
    global type_files
    global path

    for name in type_files : 
        name_files = f"{path}{study}/legend/{name}_{study}_itol.txt"
        file = open(name_files,"w")
        file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATASET_LABEL\t{name}\nDATA\n")       
        file.close()
    for name in type_files : 
        name_files = f"{path}{study}/legend/{name}_{study}NSmGn_itol.txt"
        file = open(name_files,"w")
        file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATASET_LABEL\t{name}\nDATA\n")       
        file.close()

    name_files = f"{path}{study}/legend/Label_{study}_itol.txt"
    file = open(name_files,"w")
    file.write(f"LABELS\nSEPARATOR TAB\nDATASET_LABEL\tOrganisme Label\nDATA\n")       
    file.close()

    name_files = f"{path}{study}/legend/NSmGn/Label_{study}NSmGn_itol.txt"
    file = open(name_files,"w")
    file.write(f"LABELS\nSEPARATOR TAB\nDATASET_LABEL\tOrganisme Label\nDATA\n")       
    file.close()

    name_files = f"{path}{study}/legend/Host_images{study}_itol.txt"
    file = open(name_files,"w")
    file.write(f"DATASET_IMAGE\nSEPARATOR TAB\nDATASET_LABEL\tHost with images\nCOLOR\t#000000\nSIZE_FACTOR\t0.6\nIMAGE_ROTATION\t0\nIMAGE_SHIFT_V\t0\nIMAGE_SHIFT_H\t0\nSHOW_INTERNAL\t0\nMARGIN\t400\nDATA\n")       
    file.close()

    name_files = f"{path}{study}/legend/NSmGn/Host_images{study}NSmGn_itol.txt"
    file = open(name_files,"w")
    file.write(f"DATASET_IMAGE\nSEPARATOR TAB\nDATASET_LABEL\tHost with images\nCOLOR\t#000000\nSIZE_FACTOR\t0.2\nIMAGE_ROTATION\t0\nIMAGE_SHIFT_V\t0\nIMAGE_SHIFT_H\t0\nSHOW_INTERNAL\t0\nMARGIN\t250\nDATA\n")       
    file.close()



   
### Function linux_line 
# This function concatenates all the species for each study into one file 

def linux_line(sp,study) :
    global path
    subprocess.run(['rm', f'{path}{study}/nucl_all_{sp}.fasta'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    subprocess.run(['rm', f'{path}{study}/prot_all_{sp}.fasta'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
    subprocess.run(f"cat {path}{study}/nucl/* >> {path}{study}/nucl_all_{sp}.fasta",shell=True)
    subprocess.run(f"cat {path}{study}/prot/* >> {path}{study}/prot_all_{sp}.fasta",shell=True)

### Function extract_id
# This function extracts the following information: Host, Host country, Strain, Nucleic identifier, Protein identifier for each genbank file.

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

        references=record.annotations.get('references')[0]
        ref=f"Authors : {references.authors};Title : {references.title};Journal : {references.journal};Pubmed ID : {references.pubmed_id};Medline ID : {references.medline_id} "

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
                strand = feature.location.strand          
            if "protein_id" in feature.qualifiers :
                ID_prot = feature.qualifiers['protein_id'][0]
        output_prot=""        
        output_prot=extract_prot_data(ID_prot,study)
        
        order_color=colour_choice(order)
        region_color=colour_choice(region)
        serogroup=find_serogroup(fullorganisme)
        color_sero=colour_choice(serogroup)

        if study=="CCHFV" : 
            organism=ID_prot
            ID=Modif_ID(str(ID_prot+strain))
            Description_nu = f"{ID}"
            Description_prot = f"{ID}"
            

        else : 
            ID=Modif_ID(str(organism+ID_prot))
            Description_nu = f"{ID}"
            Description_prot = f"{ID}"
        
        Description_nu.replace(" ","")
        Description_prot.replace(" ","")



        write_file(Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,sp,serogroup,color_sero)
        
        
        line=f"{fullorganisme}\t{ID}\t{taxonomy}\t{serogroup}\t{ID_prot}\t{sizeprot}\t{output_prot}{Nu}\t{sizenu}\t{strand}\t{strain}\t{region}: {localisation}\t{order}: {host}\t{ref}\n"
        line_prot_=f"{ID}\t{output_prot}"
        line_prot_.replace(">","")
        line_prot_.replace("<","")
        line_prot_.replace("GP38","")
        line_prot_.replace("Gn","")

        name=f"{path}Metadata_studies.txt"
        allmetadata=open(name,'a')
        allmetadata.write(line)
        allmetadata.close()

        name=f"{path}Prot_file_size.txt"
        protfile=open("name","a")
        protfile.write(line_prot_)
        protfile.close()

    return strain,ID_prot,Nu

###Function modif_name_prot
#The purpose of this function is to identify the protein present in the gp

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

##Function extract_prot_data
#This function extract the position of each proteins of each organism in order to save everything in the metadatafiles

def extract_prot_data(ID_prot,study) : 
    global mail
    all_proteins=["MLD","GP38","Gn","NSm","Gc"]
    Entrez.email = mail
    handle = Entrez.efetch(db="protein", id=ID_prot, rettype="gp", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    time.sleep(2)  
    with open(f"{path}{study}/prot/{ID_prot}.gp", "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")
        output_handle

    record = SeqIO.read(f"{path}{study}/prot/{ID_prot}.gp", "genbank")
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
    subprocess.run(['rm', f'{path}{study}/prot/{ID_prot}.gp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    return output

#The serogroup file 
df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0)
df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)

##Function find_serogroup 
#This function 
def find_serogroup(organism) : 
    global df
    global path
    end_line=df.shape[0]
    for line in range(0,end_line) : 
        if df.iloc[line,0].lower() in organism.lower() or organism.lower() in df.iloc[line,0].lower() or organism.replace(" ","").lower() in df.iloc[line,0].replace(" ","").lower() : 
            return str(df.iloc[line,1])
    name=f"{path}Serogroup.txt"
    file = open(name,"a")
    file.write(f"{organism}|Na\n")       
    file.close()
    df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)
    find_serogroup(organism)


#######Function for the Script mark_prot##################


##Function compt_ali_nucl
# This function uses the position of the proteins of the target organism (proteic sequence) to identify the position of the proteins in the alignment.
# It is important not to just translate the protein sequence into RNA in order to keep the intergenic regions.

def compt_ali_nucl(sequence,debseq,pos_prot) :
    pos_align=[]
    compt_ali=0
    compt_nu=0
    compt_AA=0
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
                    compt_AA+=1
                    if compt_AA in pos_prot :
                            
                        if deb==True : 
                            deb=False
                            pos_align.append(compt_ali-2)
                            pos_prot.remove(compt_AA)
                        else : 
                            deb=True
                            pos_align.append(compt_ali)
                            pos_prot.remove(compt_AA) 
                        
                    codon=""
    return pos_align
  
##Function compt_alignement
# This function uses the position of the target organism's proteins (sequence proteic) to identify the position of the proteins in the alignment.

def compt_alignement (sequence,pos_prot) : 
    
    pos_align=[]
    compt_ali=0
    compt_AA=0
    mark_prot=""
    deb=True
    for aa in sequence : 
        if aa=='-' : 
            compt_ali+=1
        else : 
            compt_ali+=1
            compt_AA+=1

        if compt_AA in pos_prot :
            if deb==True : 
                mark_prot+="|"
                mark_prot+=aa
                deb=False
                pos_align.append(compt_ali)
                pos_prot.remove(compt_AA)
            else : 
                mark_prot+=aa
                mark_prot+="|"
                pos_align.append(compt_ali)
                deb=True
                pos_prot.remove(compt_AA)
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
            name=f"{path}{study}{res_file}"+str(Protname[index])+sp+".fasta"
            file = open(name,"a")
            text=str (">"+ID+'\n'+Seq(prot)+'\n')
            file.write(text)
            file.close()
            
#######Function for the Script delete_dbl##################
##Function select_doublon
#This function identified all the possible duplicate sequences
#The possible duplicates are identified using a specific threshold; this threshold was found using the dnadist file and the visible rake of the phylogenetic tree (specific to each study and not provided by this script). 

def select_doublon(protein,study,target,ths,OG_list) :
    global path
    df = pd.read_csv(f"{path}{study}{res_file}{protein}nucl_dnadist.txt", delimiter="\t",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(float)
    
    
    Doublon=onlygap(protein,study,target,OG_list)
    end_line=df.shape[0]
    end_col=df.shape[1]
    line=0
    col=0
    while line <end_line:
        col=line+1
        while col <end_col-1 :
            if df.iloc[line,col] <=ths and df.iloc[line,col]!=-1  :
                if target not in df.index[line] and (df.index[line] not in OG_list or any(og not in df.index[line] for og in OG_list)):
                    ID=str(df.index[line])
                    Doublon.append(ID)
                    df=df.drop(columns=ID)
                    df=df.drop(index=ID)
                    end_line=df.shape[0]
                    end_col=df.shape[1]
                    line=0

            col+=1
            

        line+=1
    noms_colonnes = df.columns
    return Doublon,noms_colonnes

def onlygap(protein,study,target,OG_list):
    Doublon=[]
    types_seq=["prot","nucl"]
    for type_seq in types_seq:
        for record in SeqIO.parse(f"{path}{study}/{res_file}{protein}{type_seq}_align.aln-fasta","fasta") :
            if set(record.seq)=={'-'} and record.id != target and (record.id not in OG_list or any(og not in record.id for og in OG_list)):
                if record.id not in Doublon :
                    Doublon.append(record.id)
    return Doublon

##Function col_name
#This function delete the duplicate sequence name 
    
def col_name(old_col_names,col_name,Protein) :

    if len(old_col_names) >0 :
        if Protein=="NSmGn" : 
            col_name=[element[-9:] for element in col_name]

        old_col_names = [element for element in old_col_names if element in col_name]
    else : 
        old_col_names=col_name

    return old_col_names

##Function delete_dbl
#This function delete the identified duplicates for each protein for each study

def delete_dbl(Doublon,study,protein,old_col_names):
    global path
    
    print(protein)
    name=f"{path}{study}{res_file}{protein}nucl_sdbl.fasta"
    
    subprocess.run(f"rm {name}",shell=True)
    if study!="CCHFV" :
        file_nucl=f"{path}{study}{res_file}{protein}nucl_trimal.fasta"
        file_prot=f"{path}{study}{res_file}{protein}prot_trimal.fasta"
    else : 
        file_nucl=f"{path}{study}{res_file}{protein}nucl_align.aln-fasta"
        file_prot=f"{path}{study}{res_file}{protein}prot_align.aln-fasta"

    for record in SeqIO.parse(f"{file_nucl}","fasta") :
        
        if record.description not in Doublon and record.description in old_col_names :
            file = open(name,"a")
            text=str (">"+record.description +'\n'+record.seq+'\n')
            file.write(text)
            file.close()


    name=f"{path}{study}{res_file}{protein}prot_sdbl.fasta"
    subprocess.run(f'rm {name}',shell=True)
    for record in SeqIO.parse(f"{file_prot}","fasta") :
        
        if record.description not in Doublon and record.description in old_col_names :
            file = open(name,"a")
            text=str (">"+record.description+'\n'+record.seq+'\n')
            file.write(text)
            file.close()


    

