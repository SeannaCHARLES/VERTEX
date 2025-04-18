from Bio import SeqIO
from Bio import Entrez
import re
import subprocess

import random
import matplotlib.pyplot as plt

import requests
cmap=plt.get_cmap("hsv")
colors=[cmap(i/100)for i in range(100)]
couleur=['#%02x%02x%02x'%(int(r*255),int(g*255),int(b*255))for r,g,b,_ in colors]
random.shuffle(couleur)
metadata={}
def choix_couleur(region):
    global couleur
    if region not in metadata : 
        choix=couleur[0]
        couleur.remove(choix)
        metadata[region]=[choix]
        return(choix)
    else : 
        return(metadata[region][0])
def region_geographique(pays) : 
    url = f"https://restcountries.com/v3.1/name/{pays}"
    response = requests.get(url)
    if response.status_code==200 : 
        data=response.json()
        if isinstance(data,list) and "region" in data[0] : 
            return data[0]["region"]
    return pays

def order(species) : 
    Entrez.email="charles.seanna@outlook.fr"
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
    


name_geol="/home/seannacharles/Documents/CCHFV/Geolocalisation_CCHFV_itol.txt"      
geoloc = open(name_geol,"w")
geoloc.write(f"TREE_COLORS\nSEPARATOR TAB\nDATA\n")       
geoloc.close()

name_host="/home/seannacharles/Documents/CCHFV/Host_CCHFV_itol.txt"      
Host_f = open(name_host,"w")
Host_f.write(f"TREE_COLORS\nSEPARATOR TAB\nDATA\n")       
Host_f.close()

with open("/home/seannacharles/Documents/CCHFV/sequence_CCHFV_segM.gb","r") as Datas:
    for record in SeqIO.parse(Datas,"genbank"): 
        Nu= record.id
        Sequence_nu = record.seq
        organisme= record.annotations.get('organism',"Inconnu")
        host="_Na_"
        location="_Na_"
        souche="_Na_"
        for feature in record.features : 
            if "host" in feature.qualifiers : 
                host = feature.qualifiers['host'][0]
            if "country" in feature.qualifiers :
                location = feature.qualifiers['country'][0]
                location=location.split(":")[0]
            if "strain" in feature.qualifiers :
                souche = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['strain'][0])
            elif "isolate" in feature.qualifiers :
                souche = re.sub(r"[/\\?%*:|\"<>]","|",feature.qualifiers['isolate'][0])
            if feature.type == "CDS" and "translation" in feature.qualifiers :
                Sequence_aa=feature.qualifiers["translation"][0]
                if "protein_id" in feature.qualifiers :
                    ID_prot = feature.qualifiers['protein_id'][0]
        Description_nu = f"{Nu}-{souche}-{host}-{location}-{organisme}"
        Description_prot = f"{ID_prot}-{souche}-{host}-{location}-{organisme}"
    
        geoloc=open(name_geol,'a')
        region=region_geographique(location)
        colors=choix_couleur(region)
        geoloc.write(f"##{location}\n")
        text=str (f"{Description_nu}\trange\t{colors}\t{region}\n")
        geoloc.write(text)
        geoloc.close()

        Host_f=open(name_host,'a')
        order_=order(host)
        colors=choix_couleur(order_)
        Host_f.write(f"##{host}\n")
        text=str (f"{Description_nu}\trange\t{colors}\t{order_}\n")
        Host_f.write(text)
        Host_f.close()


