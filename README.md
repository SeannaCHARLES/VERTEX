# Script for phylogenetic studies
## Project description
Le but final de ce script est de générer des arbres phylogénétique afin que vous puissiez les analyser et les interprétés.
Concraitement ce script permet de  :  
* Download all the sequences in the Genbank format.

    This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.
* Seperate all the nucleic sequences and the proteic sequences into single files
    
    From the Genbank files the nucleic and proteic sequences are being extract and save in 2 folder M_nucl and M_prot (M that represent the segment)
* Create a **Metadata_files.txt** that contains different informations 
        
    This step download all the sequences in the gp format and delete it after
* Mark the position of the different proteins in the proteic sequences and the nucleic sequences

* Seperate the different proteins and save them 
    
    one file for each proteins.
* Identified the best model (using iqtree2)

* Create the phylogenetic tree using the best model (using iqtree2)

* 
## Instructions to use 
## Instructions to adapt the code to your needs
This scripts is currently working for the peribunyaviridae the nairoviridae the CCHFV strains and the oropouche strains if you want to use this script on other organism you need to modifiy the following points : 
### Path and folders names
#### Modify
You need to change the original path and also the folders names.
To do so all you have to do is to change the variable *path* and the variable in the list *studies* in the *execution.sh* file.
#### Add species
If you want to add another species add the name to that species in the *studies* list in the file **execution.sh**

### Mark proteins

For this step it is important to change the different position of each sequences.and to change the target ID.
You also need to change the startseq.
The startseq was found by using blast and align the proteic sequences to the nucleic sequences. 



**If you want to add another species it is recommanded to rerun the entire script in order to have the same caption for all your studied species**

### If you you are not using the serogroup 
### If you want to use these script for other species
### If you only want to create the itol legend
### If you only want to create the Metadata files 
## Pré-requis
## Installation 
## Output files

### Metadata_files.txt
This file contains the following informations(if they are available in the NCBI genbank file format and the gb file format) on all the organims you choose.

        Virus name	ID itol tree	Taxonomy	Serogroup	Protein accesion number	Size Proteique sequence	Start MLD	End MLD	Sart GP38	End GP38	Start Gn	End Gn	Start NSm	End NSm	Start Gc	End Gc	Nucleique accession number	Size Nucelique sequence	Strand	Strain	Geographical region	Host	Host source	References

