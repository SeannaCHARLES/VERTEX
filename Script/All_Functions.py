#Created by CHARLES Seanna 
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


#######Function for the Script Genbank extract##################
#This part create a set of colors to create the itol legends 


#This is the diffent type of legends
type_files=["Geolocalisation","Host","Serogroup"]

### Function Modif_ID 
# It is necessary to shorten the sequence identifier for dnadist (maximum 10 characters)
# It is important to uniformize the identifiers 
# It is important to remove all special characters that could cause problems in the future use of the software.
# It is important to remove terms that do not allow the different organisms to be differentiated (Add the terms you want in this function)



###################################################################GENBANK EXTRACT PYTHON FILE###################################################################
def Modif_ID(ID) : 
    #These changes are necessary to ensure the uniformity of the identifiers used by the various organisations.
    #Add you own modification if you feel the need
    ID=ID.lower()    
    ID=ID.replace("_","").replace("-","").replace("'","").replace(".","").replace(" ","")
    ID=ID[:9]
    return ID

### Function Colour_choice
# A colour is assigned to each parameter that may appear in a legend
# The colours are chosen at random 
# All families will have the same captions

#Legencolors save for each legend label the associated colors to have the same color for the same label 
ratp_colors = [
   
    "#FFCD00", 
    "#003CA6", 
    "#837902", 
    "#6EC4E8", 
    "#CF009E", 
    "#FF7E2E", 
    "#6ECA97", 
    "#FA9ABA", 
    "#E19BDF", 
    "#B6BD00", 
    "#C9910D", 
    "#704B1C", 
    "#007852", 
    "#62259D", 

    "#E2231A", 
    "#003DA5", 
    "#F7A600", 
    "#008641", 

    "#BBBD00", 
    "#006AB6", 
    "#00A7E1", 
    "#BB77B7", 
    "#FF5A00", 
    "#6E4D8C", 
    "#DA1A32"  
]

legencolors={}
legencolors["Na"]=["#979a9a"]




def colour_choice (caption,colors_list,type_,path,study):
    global ratp_colors

    if len(colors_list)<=0 : 
        colors_list=ratp_colors.copy()

    if caption not in legencolors.keys() : 
        choice=colors_list[0]
        colors_list.remove(choice)
        legencolors[caption]=[choice]

        name_legend=f"{path}/{study}/All_legend_caption.txt"
        All_legend=open(name_legend,'a')
        text=str (f"{caption}\t{choice}\n")
        All_legend.write(text)
        All_legend.close()
        return(choice,colors_list)
    else : 

        return(legencolors[caption][0],colors_list)

### Function modif_localisation
# For each geographical region, it is important to keep only those elements that identify the host's country of origin.

def modif_localisation(localisation) :
    localisation=localisation.split(":")[0]
    return localisation

### Function geographic_region
# For the captions, each country is replaced by its geographical region.
# At the moment I can't get pycountry_convert to work, which would simplify this function.
# If the geographical region can't be find the country will be kept

def geographic_region(pays,email) : 
    latitude=0
    longitude=0

    geolocator = Nominatim(user_agent=f"{pays}continent{email}")
    location = geolocator.geocode(pays,timeout=10)

    latitude=location.latitude
    longitude=location.longitude
    if latitude==0 or longitude==0 : 
        return pays
    elif -35 <= latitude <= 37 and -20 <= longitude <= 55:
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

def modif_host(host,list_speci_order) : 

    for sp_or in list_speci_order : 
        speci_order=sp_or.split("=")
        if speci_order[0].lower() in host.lower() : 
            host=speci_order[1]

    host=host.split(";")[0]
    host=host.replace("sp.",'')
    host=host.replace("sentinel.",'')
    host=re.sub(r"\([^)]*\)", "", host) #replace all the text betwenn ()

    return host

### Function deforder 
# For the captions, each species is replaced by its order
# If the order can't be find the species will be kept
# This function send a request to the NCBI database
def extract_id_tax(record,mail) : 
    Entrez.email=mail

    tax_id = record["IdList"][0]
    tax_info= Entrez.efetch(db="taxonomy",id=tax_id,retmode="xml")
    data=Entrez.read(tax_info)[0]

    if data["Rank"]== "order" : 
        return data["ScientificName"]

    else :    
        for tax_id_group in data["LineageEx"] : 
            if tax_id_group["Rank"]=="order":
                return tax_id_group["ScientificName"]



def deforder(species,mail) : 

    Entrez.email=mail
    sp =species.split(" ")

    try : 
        search= Entrez.esearch(db="taxonomy",term=species)
        record=Entrez.read(search)
        
        if not record["IdList"]: # Pas d'indentifiant taxonomy
            
            for element in sp : 
                try : 
                    search= Entrez.esearch(db="taxonomy",term=element)
                    record=Entrez.read(search)
                
                    if record["IdList"] :

                        return extract_id_tax(record,mail)
                except Exception as e :
                    print()
            print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add it in the corresponding line of the information_file.txt such as.\n ;{species.lower()}='Scientificname'")
      
            return species

        else : 

            return extract_id_tax(record,mail)

    except Exception as e :

        for element in sp : 
            search= Entrez.esearch(db="taxonomy",term=element)
            record=Entrez.read(search)

            if record["IdList"] :
                return extract_id_tax(record,mail)

        print(f"\nThe species '{species}' was not identified by the NCBI taxonomy database\nPlease check the scientific name of this species and add it in the corresponding line of the information_file.txt such as.\n ;{species.lower()}='Scientificname'")
        return species



### Function write_file
# This function writes the caption files required for itol
# This function writes individual files for each nucleotide and protein sequence 

def write_file(path,Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,serogroup,color_sero,sp,ccprot) : 
    global type_files
    namespeprot={}
    
    for pos,element in enumerate(ccprot) : 
        if element !="" : 
            if "=" in element : 
                if len(element.split("=")[1].split("+"))>=1:
                    namespeprot[element.split("=")[0]]=element.split("=")[1].split("+")

            
        
    
    
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
    text=str (f"{ID}\trange\t{color_sero}\t{serogroup}\n{ID}\tbranch\t{color_sero}\tdashed\t7\n")
    sero.write(text)
    sero.close()

    HostI_f=open(name_host_images,'a')
    HostI_f.write(f"##{host}\n")
    text=str (f"{ID}\t-1\t1\t0\t0\t0\thttps://raw.githubusercontent.com/SeannaCHARLES/Images/main/{order.lower()}.png\n")
    HostI_f.write(text)
    HostI_f.close()
    
    if namespeprot != {}: 
        
        for speprot,values in namespeprot.items() : 
            name_geol=f"{path}{study}/legend/{type_files[0]}_{study}{speprot.lower()}_itol.txt"
            name_host=f"{path}{study}/legend/{type_files[1]}_{study}{speprot.lower()}_itol.txt"
            name_label=f"{path}{study}/legend/Label_{study}{speprot.lower()}_itol.txt"
            name_sero=f"{path}{study}/legend/{type_files[2]}_{study}{speprot.lower()}_itol.txt"
            name_host_images=f"{path}{study}/legend/Host_images{study}{speprot.lower()}_itol.txt"
            name_range=f"{path}{study}/legend/Range_{study}{speprot.lower()}_itol.txt"


            geoloc=open(name_geol,'a')
            geoloc.write(f"\n##{localisation}\n")
            text="\n".join(f"{element}{ID}\trange\t{region_color}\t{region}" for element in values)
            geoloc.write(text)
            geoloc.close()

            Host_f=open(name_host,'a')
            Host_f.write(f"\n##{host}\n")
            text="\n".join(f"{element}{ID}\trange\t{order_color}\t{order}" for element in values)
            Host_f.write(text)
            Host_f.close()

            label=open(name_label,'a')
            label.write(f"\n##{organism}\n")
            text="\n".join(f"{element}{ID}\t{element}_{organism}" for element in values)
            label.write(text)
            label.close() 

            sero=open(name_sero,'a')
            sero.write(f"\n##{serogroup}\n")
            text="\n".join(f"{element}{ID}\trange\t{color_sero}\t{serogroup}\n{element}{ID}\tbranch\t{color_sero}\tdashed\t7" for element in values)
            sero.write(text)
            sero.close()
    
            label=open(name_range,'a')
            label.write(f"\n##{organism}\n")
            text="\n".join(f"{element}{ID}\t{color_sero}" for element in values)
            label.write(text)
            label.close() 

            HostI_f=open(name_host_images,'a')
            HostI_f.write(f"\n##{host}\n")
            text="\n".join(f"{element}{ID}\t-1\t1\t0\t0\t0\thttps://raw.githubusercontent.com/SeannaCHARLES/Images/main/{order.lower()}.png"for element in values)
            HostI_f.write(text)
            HostI_f.close()

    nomprot_file=f"{Nu}_prot"
    nomprot_file=nomprot_file.replace(" ","")
    nomnucl_file=f"{Nu}_nucl"
    nomnucl_file=nomnucl_file.replace(" ","")


    if sp =="OGspe" :         
        with open (f"{path}{study}/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")  
        with open (f"{path}{study}/nucl/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/prot/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")
    else : 
        with open (f"{path}{study}/nucl/{nomnucl_file}.fasta","w") as nu : 
            nu.write(f">{ID}\n{Sequence_nu}\n")
        with open (f"{path}{study}/prot/{nomprot_file}.fasta","w") as aa : 
            aa.write(f">{ID}\n{Sequence_aa}\n")

  

### Function first_call 
# This function writes the header for caption files for Itol
      
def first_call(path,study,cop) : 
    global type_files
    subprocess.run(['rm', f'{path}{study}/legend/*'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    namespeprot=cop.copy()
    for pos,element in enumerate(namespeprot) : 
        namespeprot[pos]=element.split("=")[0]
 
    for name in type_files : 
        name_files = f"{path}{study}/legend/{name}_{study}_itol.txt"
        file = open(name_files,"w")
        file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATASET_LABEL\t{name}\nDATA\n")       
        file.close()
    for name in type_files : 
        for speprot in namespeprot : 
            name_files = f"{path}{study}/legend/{name}_{study}{speprot.lower()}_itol.txt"
            file = open(name_files,"w")
            file.write(f"TREE_COLORS\nSEPARATOR TAB\nDATASET_LABEL\t{name}\nDATA\n")       
            file.close()

    name_files = f"{path}{study}/legend/Label_{study}_itol.txt"
    file = open(name_files,"w")
    file.write(f"LABELS\nSEPARATOR TAB\nDATASET_LABEL\tOrganisme Label\nDATA\n")       
    file.close()

    name_files = f"{path}{study}/legend/Host_images{study}_itol.txt"
    file = open(name_files,"w")
    file.write(f"DATASET_IMAGE\nSEPARATOR TAB\nDATASET_LABEL\tHost with images\nCOLOR\t#000000\nSIZE_FACTOR\t0.6\nIMAGE_ROTATION\t0\nIMAGE_SHIFT_V\t0\nIMAGE_SHIFT_H\t0\nSHOW_INTERNAL\t0\nMARGIN\t400\nDATA\n")       
    file.close()

    if namespeprot != {}: 
        for speprot in namespeprot : 
            name_files = f"{path}{study}/legend/Label_{study}{speprot.lower()}_itol.txt"
            file = open(name_files,"w")
            file.write(f"LABELS\nSEPARATOR TAB\nDATASET_LABEL\tOrganisme Label\nDATA\n")       
            file.close()

            name_files = f"{path}{study}/legend/Host_images{study}{speprot.lower()}_itol.txt"
            file = open(name_files,"w")
            file.write(f"DATASET_IMAGE\nSEPARATOR TAB\nDATASET_LABEL\tHost with images\nCOLOR\t#000000\nSIZE_FACTOR\t0.2\nIMAGE_ROTATION\t0\nIMAGE_SHIFT_V\t0\nIMAGE_SHIFT_H\t0\nSHOW_INTERNAL\t0\nMARGIN\t250\nDATA\n")       
            file.close()

            name_files = f"{path}{study}/legend/Range_{study}{speprot.lower()}_itol.txt"
            file = open(name_files,"w")
            file.write(f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL	RectZoneExample\tCOLOR\t#000000\nDATA\n")       
            file.close()




   
### Function concatenates 
# This function concatenates all the species for each study into one file 

def concatenates(path,study) :

    subprocess.run(['rm', f'{path}{study}/nucl_all_{study}.fasta'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    subprocess.run(['rm', f'{path}{study}/prot_all_{study}.fasta'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
    subprocess.run(f"cat {path}{study}/nucl/*.fasta >> {path}{study}/nucl_all_{study}.fasta",shell=True)
    subprocess.run(f"cat {path}{study}/prot/*.fasta >> {path}{study}/prot_all_{study}.fasta",shell=True)
    subprocess.run(f"cat {path}{study}/nucl/*.gp >> {path}{study}/nucl_all_{study}.gp",shell=True)
    subprocess.run(f"cat {path}{study}/prot/*.gp >> {path}{study}/prot_all_{study}.gp",shell=True)
### Function extract_id
# This function extracts the following information: Host, Host country, Strain, Nucleic identifier, Protein identifier for each genbank file.
region_colors=ratp_colors.copy()
serogroup_colors=ratp_colors.copy()
host_colors=ratp_colors.copy()

def extract_id(path,Datas,study,sp,list_prot,email,Diff_name,metadata_name_prot,speci_order,dicy,namepeprot) :
    global region_colors
    global serogroup_colors
    global host_colors
    global allmetadata


    namespeprot=namepeprot.copy()
 

    for record in SeqIO.parse(Datas,"genbank"): 
        #Variable initialisation 
        host="Na"
        localisation="Na"
        strain="Na" 
        order="Na"
        region="Na"
        ID=''
        serogroup='Na'
        sizeprot=0
        Nu= record.id
        Sequence_nu = record.seq
        sizenu=len(Sequence_nu)
        organism= record.annotations.get('organism',"Inconnu").replace("virus","")
        fullorganisme=organism

        Sequence_aa=""
        ID_prot=[]
        taxonomy=record.annotations.get('taxonomy')[7:]
        taxonomy=' : '.join(taxonomy)

        references=record.annotations.get('references')[0]
        ref=f"Authors : {references.authors};Title : {references.title};Journal : {references.journal};Pubmed ID : {references.pubmed_id};Medline ID : {references.medline_id} "

        for feature in record.features : 
            if "host" in feature.qualifiers : 
                host = feature.qualifiers['host'][0]
                host=modif_host(host,speci_order)
                order = deforder(host,email) 

            if "country" in feature.qualifiers :
                localisation = feature.qualifiers['country'][0]
                localisation=modif_localisation(localisation)
                region = geographic_region(localisation,email)

            elif "geo_loc_name" in feature.qualifiers :
                localisation = feature.qualifiers['geo_loc_name'][0]
                localisation=modif_localisation(localisation)
                region = geographic_region(localisation,email)

            if "strain" in feature.qualifiers :
                strain = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['strain'][0])

            elif "isolate" in feature.qualifiers :
                strain = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['isolate'][0])

            if feature.type == "CDS" and "translation" in feature.qualifiers : 
                strand = feature.location.strand          
            if "protein_id" in feature.qualifiers :
                ID_prot.append(feature.qualifiers['protein_id'][0])
        

        if not ID_prot : 
            print(f"This ID : {Nu} does not have any proteins ID It will not be inclued in this study")
        else :       
            output_prot,output_pos,sizeprot,Sequence_aa=extract_prot_data(ID_prot,study,list_prot,email,Diff_name,metadata_name_prot,dicy,path)
        
            order_color,host_colors=colour_choice(order,host_colors,"order",path,study)
            region_color,region_colors=colour_choice(region,region_colors,"region",path,study)
            serogroup=find_serogroup(fullorganisme,path)
            color_sero,serogroup_colors=colour_choice(serogroup,serogroup_colors,"serogroup",path,study)

            ID=Modif_ID(str(organism[:3]+Nu[3:]))
            Description_nu = f"{ID}"
            Description_prot = f"{ID}"
            
            if strand =="-" or strand =="-1" : 
                Sequence_nu=Sequence_nu[::-1]
                
            write_file(path,Description_nu,Description_prot,Sequence_aa,Sequence_nu,localisation,region,region_color,order,order_color,host,strain,Nu,ID_prot,study,organism,ID,serogroup,color_sero,sp,namespeprot)
        
        
            line=f"{study}\t{fullorganisme}\t{ID}\t{taxonomy}\t{serogroup}\t{"-".join(ID_prot)}\t{sizeprot}\t{output_prot}{Nu}\t{sizenu}\t{strand}\t{strain}\t{region}: {localisation}\t{order}: {host}\t{ref}\n"
            line_prot_=f"{study}\t{ID}\t{output_pos}\n"

            line_prot_=line_prot_.replace(":","")
            name=f"{path}{study}/Metadata_{study}.txt"
            allmetadata=open(name,'a')
            allmetadata.write(line)
            allmetadata.close()
            line_prot_=line_prot_.replace(">","")
            line_prot_=line_prot_.replace("<","")

            name=f"{path}{study}/Pos_prot_{study}.txt"
            protfile=open(name,"a")
            protfile.write(line_prot_)
            protfile.close()

    

    return strain,ID_prot,Nu

###Function modif_name_prot
#The purpose of this function is to identify the protein present in the gp

def modif_name_prot(prot,Diff_name,metadata_name_prot) :


    for prot_name in metadata_name_prot : 
        if prot_name.lower() in prot.lower() : 
            prot=prot_name.lower()

    for diff_prot in Diff_name :
        
        DP=diff_prot.split("=")
        prot=prot.lower().replace("-","").replace("\xa0", " ").replace("PROTEIN","").replace("_","").replace(" ","")
        prot=prot.lower().replace(DP[0].lower(),DP[1].lower())
        

    for prot_name in metadata_name_prot : 
        if prot_name.lower() in prot.lower() : 
            prot=prot_name.lower()
    
    return prot

##Function extract_prot_data
#This function extract the position of each proteins of each organism in order to save everything in the metadatafiles

def extract_prot_data(ID_prot_list,study,all_proteins,mail,Diff_name,metadata_name_prot,dicy,path) : 

    Entrez.email = mail
    output=""  
    output_pos=""
    proteins={}
    All_end=0
    Sequence_aa=''
    prot_names=metadata_name_prot.copy()
    
    for ID_prot in ID_prot_list : 
        
        try :
            handle = Entrez.efetch(db="protein", id=ID_prot, rettype="gp", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            time.sleep(1)  
            with open(f"{path}{study}/prot/{ID_prot}.gp", "w") as output_handle:
                SeqIO.write(record, output_handle, "genbank")
                output_handle

            record = SeqIO.read(f"{path}{study}/prot/{ID_prot}.gp", "genbank")
            
            for feature in record.features :
            
                if 'product' in feature.qualifiers : 
                    prot=(feature.qualifiers['product'][0])
                    prot=modif_name_prot(prot,Diff_name,metadata_name_prot)
                    proteins=extract_start_end(prot,proteins,feature,All_end)
                    
                    
                if 'region_name' in feature.qualifiers : 
                    prot=(feature.qualifiers['region_name'][0])
                    prot=modif_name_prot(prot,Diff_name,metadata_name_prot)
                    proteins=extract_start_end(prot,proteins,feature,All_end)
                    
                    
            Sequence_aa+=record.seq
            All_end+=len(record.seq)

            # subprocess.run(['rm', f'{path}{study}/prot/{ID_prot}.gp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        except Exception as e :
            print(f"This ID : {ID_prot} were not indentified by the NCBI database or error {e}")

    
    for prot in metadata_name_prot : 
        
        if prot.lower() in proteins  :

            output+=f'{str(proteins[prot.lower()][0])}\t{str(proteins[prot.lower()][1])}\t'
            output_pos+=f'{str(proteins[prot.lower()][0])}\t{str(proteins[prot.lower()][1])}\t'
        elif any(prot.lower() in ID_p.lower() for ID_p in metadata_name_prot): 
                
            output+=f"nan\tnan\t"
            output_pos+=f"nan\tnan\t"

    return output,output_pos,All_end,Sequence_aa

#The serogroup file 
def extract_start_end(prot,proteins,feature,All_end):
    start=(feature.location.start)+All_end+1
    end=(feature.location.end)+All_end+1
     
    if start==0 : 
        start=1
                
    if "+" in prot : 
        two_prot=prot.split("+")
                    
        if any(two_prot[1].lower() in p.lower() for p in proteins) : 
            proteins[two_prot[1]][1]=end
        else : 
            proteins[two_prot[0]]=[start,"nan"]
            proteins[two_prot[1]]=["nan",end]
                          
    else : 
        if prot.lower() in proteins : 
       
            if proteins[prot][1]< end :
                proteins[prot][1]=end

        else : 
            proteins[prot]=[start,end] 
 
    return proteins

##Function find_serogroup 
#This function 
def find_serogroup(organism,path) : 

    df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)
    global collaps

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
    return "Na"
    

    


###################################################################Function for the Script delete_dbl###################################################################

##Function select_doublon
#This function identified all the possible duplicate sequences
#The possible duplicates are identified using a specific threshold; this threshold was found using the dnadist file and the visible rake of the phylogenetic tree (specific to each study and not provided by this script). 

def select_doublon(path,study,ths,keep_list,target) :

    df = pd.read_csv(f"{path}{study}/{study}_dnadist.txt", delimiter="\t",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(float)


    Doublon=onlygap(path,study,target,keep_list)
    end_line=df.shape[0]
    end_col=df.shape[1]
    line=0
    col=0
    while line <end_line:
        col=line+1
        while col <end_col-1 :
            if df.iloc[line,col] <=ths and df.iloc[line,col]!=-1  :
                if (df.index[line] not in target[9:].lower() or target[9:].lower() not in df.index[line])  and (df.index[line] not in keep_list or any(prot[9:].lower() not in df.index[line] for prot in keep_list)): #or any(og not in df.index[line] for og in keep_list):
                    ID=str(df.index[line])
                    Doublon.append(ID)
                    df=df.drop(columns=ID)
                    df=df.drop(index=ID)
                    end_line=df.shape[0]
                    end_col=df.shape[1]
                    line=0

            col+=1
            

        line+=1
    return Doublon

##Function delete_dbl
#This function delete the identified duplicates for each protein for each study

def delete_dbl(path,Doublon,study):
    types_seq=["nucl",'prot']
    for type_seq in types_seq :
    
        name=f"{path}{study}/{type_seq}_sdbl_{study}.fasta"
    
        subprocess.run(f"rm {name}",shell=True)

        for record in SeqIO.parse(f"{path}{study}/{type_seq}_all_align_{study}.aln-fasta","fasta") :
        
            if record.id not in Doublon :
                file = open(name,"a")
                text=str (">"+record.id +'\n'+record.seq+'\n')
                file.write(text)
                file.close()


def onlygap(path,study,target,keep_list):
    Doublon=[]
    types_seq=["prot","nucl"]
    for type_seq in types_seq:
        for record in SeqIO.parse(f"{path}{study}/nucl_all_align_{study}.aln-fasta","fasta") :
            if set(record.seq)=={'-'} and (record.id not in target[9:].lower() or target[9:].lower() not in record.id) and (record.id not in keep_list or any(prot[9:].lower() not in record.id for prot in keep_list)):
                if record.id not in Doublon :
                    Doublon.append(record.id)
    return Doublon


###################################################################Function for the Script Mark_script###################################################################

##Function compt_ali_nucl
# This function uses the position of the proteins of the target organism (proteic sequence) to identify the position of the proteins in the alignment.
# It is important not to just translate the protein sequence into RNA in order to keep the intergenic regions.

def nucleic_sequences(study,ID,sequence,debseq,proteins,sp,resfile,path) :

    Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
    ligne = Pos_prot_final[Pos_prot_final.iloc[:, 0] == ID]
    position_proteines=[]
    end_col=Pos_prot_final.shape[1]
    for col in range (1,end_col-1,2) : #Pos_prot_final.columns:
        start = ligne.iloc[0,col]
        end = ligne.iloc[0,col+1]
        position_proteines.append([start,end])
    #print(position_proteines)
    prot_position=dict(zip(proteins, position_proteines))
    for prot,pos in prot_position.items() : 
        #print(f"nucleic seq & OG : {prot}:{pos} : {not any("nan" in str(element) for element in pos)} ! {not any("ABS" in str(element) for element in pos)}!{pos!="nan"}")
        if not any("ABS" in str(element) for element in pos) and not any("nan" in str(element) for element in pos) : 
            pos = [int(float(x)) for x in pos]
            name=f"{path}{study}{resfile}{prot.lower()}_{sp}_nucl.fasta"
            mark_prot=""
            compt_AA=0
            compt_nu=0
            compt_codon=0
            inprot=False
    
            for nu in sequence : 

                if inprot==True : 
                    mark_prot+=nu
                if nu!='-' : 
                    compt_nu+=1
                    if compt_nu >= debseq :
                        compt_codon+=1
                
                        if compt_codon==3 : 
                            compt_codon=0
                            compt_AA+=1

                            if compt_AA in pos :
                                mark_prot+=nu
                                pos.remove(compt_AA)
                                if inprot==False : 
                                    inprot=True
                                else : 
                                    inprot=False
            
            file = open(f"{name}","a")
            file.write(f">{ID}\n{mark_prot}\n")
            file.close()
 
def nucleic_sequences_OG(study,ID,sequence,debseq,pos_prot,proteins,sp,resfile,path) :

    position_proteines=[pos_prot[i:i+2] for i in range(0, len(pos_prot), 2)]
    prot_position=dict(zip(proteins, position_proteines))

    for prot,pos in prot_position.items() : 
        #print(f"nucleic seq & OG : {prot} :{pos}: {not any("nan" in str(element) for element in pos)} ! {not any("ABS" in str(element) for element in pos)}!{pos!="nan"}")
        if not any("ABS" in str(element) for element in pos) and not any("nan" in str(element) for element in pos) : 
            pos = [int(float(x)) for x in pos]
            name=f"{path}{study}{resfile}{prot.lower()}_{sp}_nucl.fasta"
            mark_prot=""
            compt_AA=0
            compt_nu=0
            compt_codon=0
            inprot=False
    
            for nu in sequence : 

                if inprot==True : 
                    mark_prot+=nu
                if nu!='-' : 
                    compt_nu+=1
                    if compt_nu >= debseq :
                        compt_codon+=1
                
                        if compt_codon==3 : 
                            compt_codon=0
                            compt_AA+=1

                            if compt_AA in pos :
                                mark_prot+=nu
                                pos.remove(compt_AA)
                                if inprot==False : 
                                    inprot=True
                                else : 
                                    inprot=False
            
            file = open(f"{name}","a")
            file.write(f">{ID}\n{mark_prot}\n")
            file.close()
##Function compt_alignement
# This function uses the position of the target organism's proteins (sequence proteic) to identify the position of the proteins in the alignment.

def proteic_sequence (study,ID,sequence,proteins,sp,resfile,path) : 
 
    Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
    ligne = Pos_prot_final[Pos_prot_final.iloc[:, 0] == ID]
    position_proteines=[]
    end_col=Pos_prot_final.shape[1]
    for col in range (1,end_col-1,2) : #Pos_prot_final.columns:
        start = ligne.iloc[0,col]
        end = ligne.iloc[0,col+1]
        position_proteines.append([start,end])
     
 
    #position_proteines=[pos_prot[i:i+2] for i in range(0, len(pos_prot), 2)]
    prot_position=dict(zip(proteins, position_proteines))
 


    for prot,pos in prot_position.items() : 
        
        name=f"{path}{study}{resfile}{prot.lower()}_{sp}_prot.fasta"
        if not any("ABS" in str(element) for element in pos) and not any("nan" in str(element) for element in pos): 
            pos = [int(float(x)) for x in pos]
            mark_prot=""
            compt_AA=0
            inprot=False
            for aa in sequence : 
                if aa!='-' : 
                    compt_AA+=1

                if inprot==True : 
                    mark_prot+=aa

                if compt_AA in pos :
                   
                    mark_prot+=aa
                    pos.remove(compt_AA)
                    if inprot==False : 
                        inprot=True
                    else : 
                        inprot=False
            file = open(f"{name}","a")
            file.write(f">{ID}\n{mark_prot}\n")
            file.close()



#Function mark_sequence
#This function marks the start and end of proteins in protein and nucleic acid sequences. 

def final_pos_IO(study,prot,study_prot,proteins,concern_study,path) : 
    
   
    chain_prot=""
    df = pd.read_csv(f"{path}{study}/Pos_prot_{study}.txt", delimiter="\t",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)
    end_line=df.shape[0]
    end_col=df.shape[1]
    for line in range (1,end_line-1,1) :
        for col in range(2,end_col-1,1) : 
  
            if df.iloc[line,col]=="nan" : 
                

                break
                #subprocess.run(f"blastp  -db -query -out -max_target_seqs -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"")


    return (df.iloc[line,col],df.iloc[line,col])



def final_pos(study,prot,study_prot,proteins,path) : 
    chain_prot=""
    df = pd.read_csv(f"{path}{study}/Pos_prot_{study}.txt", delimiter="\t",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)
    end_line=df.shape[0]
    end_col=df.shape[1]

    pos_col_prot={}
    for col in df.columns : 
        for item in prot :
            item=item.split(";")
            for pos_prot in item :
                if pos_prot not in pos_col_prot :
                    if pos_prot in col :
                        pos_col_prot[pos_prot]=[df.columns.get_loc(col),df.columns.get_loc(col)+1]
    final_pos_prot={}
   
    for line in range (0,end_line) :
  
        pos_col=col_pos(str(df.index[line]),pos_col_prot,study_prot)
        final_pos_prot[df.iloc[line,0]]=[] 
        for col, prot in zip(range(1,end_col-1), proteins[str(df.index[line])]):

            if df.iloc[line,col]=="nan" :
                if col%2!=0   : 
                    if col >1 :
                        if df.iloc[line,col-1]!="nan":
                            final_pos_prot[df.iloc[line,0]].append(int(float(df.iloc[line,col-1]))+1) 
                        else :
                            final_pos_prot[df.iloc[line,0]].append(int(prot)) 
                    else : 
                        final_pos_prot[df.iloc[line,0]].append(int(prot)) 
                else : 
                    if col+1<end_col :
                        if df.iloc[line,col+1]!="nan":
                            final_pos_prot[df.iloc[line,0]].append(int(float(df.iloc[line,col+1]))-1) 
                        else : 
                            final_pos_prot[df.iloc[line,0]].append(int(prot))
                    else : 
                        final_pos_prot[df.iloc[line,0]].append(int(prot)) 

            else :
                final_pos_prot[df.iloc[line,0]].append(int(float(df.iloc[line,col])))
   

        
        for c in range (1,len(final_pos_prot[df.iloc[line,0]] )-1) :
            if c%2 ==0 : 
                if final_pos_prot[df.iloc[line,0]][c+1] <= final_pos_prot[df.iloc[line,0]][c] : 
                    final_pos_prot[df.iloc[line,0]][c]=final_pos_prot[df.iloc[line,0]][c+1]-1
                if final_pos_prot[df.iloc[line,0]][c-1] >= final_pos_prot[df.iloc[line,0]][c] : 
                    final_pos_prot[df.iloc[line,0]][c]=final_pos_prot[df.iloc[line,0]][c-1]+1
            else : 
                if final_pos_prot[df.iloc[line,0]][c-1] >= final_pos_prot[df.iloc[line,0]][c] : 
                    final_pos_prot[df.iloc[line,0]][c]=final_pos_prot[df.iloc[line,0]][c-1]+1
                if final_pos_prot[df.iloc[line,0]][c+1] <= final_pos_prot[df.iloc[line,0]][c] : 
                    final_pos_prot[df.iloc[line,0]][c]=final_pos_prot[df.iloc[line,0]][c+1]-1

    ligne_chaine = '\t'.join(df.iloc[[0]].astype(str))
    name=f"{path}{study}/Pos_prot_final_{study}.fasta"
    file = open(f"{name}","w")
    file.write(f"{ligne_chaine}\n")
    file.close()
    return final_pos_prot,name

def col_pos(nm_study,pos_col_prot,study_prot) :

    range_col=[]
    for element in study_prot[nm_study]:
        for e in pos_col_prot[element]:
            range_col.append(int(e))

    return(range_col)



###################################################################Function for the Script WIWO###################################################################


def deleteIN(study,IO,protein,organismes,res_file,types_seq,path) : 

    for type_seq in types_seq : 
        subprocess.run(f" cp {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta {path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta" ,shell=True)
        subprocess.run(["rm", f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"] ,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,shell=True)
        filename=f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"

        for record in SeqIO.parse(f"{path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta","fasta") :
            if any(org[:6] in record.id for org in organismes):
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()

def deleteOUT(study,IO,protein,organismes,res_file,types_seq,path) : 

    for type_seq in types_seq : 
        subprocess.run(f"cp {path}{study}{res_file}{protein}_{study}_{type_seq}.fasta {path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta",shell=True)
        subprocess.run(["rm", f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"], stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,shell=True)
        filename=f"{path}{study}{res_file}{protein}_{study}_{type_seq}.fasta"

        for record in SeqIO.parse(f"{path}{study}{res_file}{protein}_{study}_{type_seq}_bef_WIWO.fasta","fasta") :
            if not any(org[:6] in record.id for org in organismes): 
                file = open(filename,"a")
                file.write(f">{record.id}\n{record.seq}\n")       
                file.close()

    


##Function Extract_data
# This function separates the different proteins that have been tagged for all sequences.

def Extract_data(file,sp,Protname,study,path) : 
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
    write_file_prot(Data,sp,Protname,study,path)


## Function write_file_prot
# This function writes a file for each protein with all the sequences studied.  
 
def write_file_prot(Data,sp,Protname,study,path) : 
 
    for ID,seq in Data.items() : 
        for index,prot in enumerate(seq) : 
            name=f"{path}{study}{res_file}"+str(Protname[index])+sp+".fasta"
            file = open(name,"a")
            text=str (">"+ID+'\n'+Seq(prot)+'\n')
            file.write(text)
            file.close()
            

def Blastp_run(DB,File,output) :
    subprocess.run(f"rm -f {output}.txt",shell=True)
    subprocess.run(f"blastp -db {DB}_DB -query {File} -outfmt 6 -out {output}.txt",shell=True)

    
def add_new_pos(protein,output) :
    prot=protein.replace("Start","").replace("End","")
    if "Start" in protein :
        missing_pos=7
    else : 
        missing_pos=8


    cmd=f"awk '/{prot}/ && (NR==1 || $4 > max) {{ max = $4; line = ${missing_pos} }} END {{ print line }}' {output}.txt"
   
    pos=subprocess.run(f"{cmd}",capture_output=True,shell=True, text=True)
   
    if pos.stdout=="\n" or pos.stdout=='' : 
        return int(-1)
    else  :
        return int(pos.stdout)

def change_pos(line,end_col,Pos_prot_final,path,study) : 
    for col in range (1,end_col) :
        if "Start" in str(Pos_prot_final.columns[col]) : 
           
            if "str" not in str(type(Pos_prot_final.iloc[line,col]))   :

                if Pos_prot_final.iloc[line,col]==-1 :
                   
                    Pos_prot_final.iloc[line,col]="nan"
                    Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and  Pos_prot_final.iloc[line,col+1]!=-1:
                    if Pos_prot_final.iloc[line,col] > Pos_prot_final.iloc[line,col+1]  : 
                        if col ==1 : 
                            Pos_prot_final.iloc[line,col]=1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                           
                        else :
                           
                            Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col-1])+1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
            elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and "str" not in str(type(Pos_prot_final.iloc[line,col-1])) : 
                Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col-1])+1
                Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                        
        else : 
            if "str" not in str(type(Pos_prot_final.iloc[line,col]))  :

                if Pos_prot_final.iloc[line,col]==-1 :
             
                    Pos_prot_final.iloc[line,col]="nan"
                    Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                elif "str" not in str(type(Pos_prot_final.iloc[line,col-1])) and  Pos_prot_final.iloc[line,col-1]!=-1:
                    if  Pos_prot_final.iloc[line,col-1] > Pos_prot_final.iloc[line,col] :
                        if col !=end_col : 
                            Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col+1])-1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                        else : 
                            Pos_prot_final.iloc[line,col]="END"
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
            elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and "str" not in str(type(Pos_prot_final.iloc[line,col-1])) : 
                Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col+1])-1
                Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
  
    
