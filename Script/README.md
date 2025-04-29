# Script for phylogenetic studies
## Project description
Le but final de ce script est de générer des arbres phylogénétique afin que vous puissiez les analyser et les interprétés.
Concraitement ce script permet de  :  
* Download all the sequences in the Genbank format.

    This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.
* Seperate all the nucleic sequences and the proteic sequences into single files
    
    From the Genbank files the nucleic and proteic sequences are being extract and save in 2 folder M_nucl and M_prot (M that represent the segment)
* Create a Metadata_files.txt that contains different informations such as  ID nucleic, ID proteic , taxonomy, and the position of the differents proteins.
        
    This step download all the sequences in the gp format and delete it after
* Mark the position of the different proteins in the proteic sequences and the nucleic sequences

    For this step it is important to change the different position of each sequences and to change the target ID.
    You also need to change the startseq 
    
    The startseq was found by using blast and align the proteic sequences to the nucleic sequences 

* Seperate the different proteins and save them 
    
    one file for each proteins.
* Identified the best model (using iqtree2)

* Create the phylogenetic tree using the best model 

* 
## Instructions to use 
## Instructions to adapt the code to your needs
This scripts is currently working for the peribunyaviridae the nairoviridae and CCHFV strains if you want to use this script on other organism you need to modifiy the following points
### If you you are not using the serogroup 
### If you want to use these script for other species
### If you only want to create the itol legend
### If you only want to create the Metadata files 
## Pré-requis
## Installation 
## Output files

### Metadata_files.txt
This file contains the following informations(if they are available in the NCBI genbank file format and the gb file format) on all the 
From these files multiple informations such as the taxonomy, , ect are being saved in a Metadata_files.txt





Output files 

