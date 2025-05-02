# Script for phylogenetic studies
## Project description
Le but final de ce script est de générer des arbres phylogénétique afin que vous puissiez les analyser et les interprétés.
Concraitement ce script permet de  :  
* Download all the sequences in the Genbank format.

    This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.

* Seperate all the nucleic sequences and the proteic sequences into single files
    
    From the Genbank files the nucleic and proteic sequences are being extract and save in 2 folder nucl and prot (M that represent the segment)

* Create a **Metadata_files.txt** that contains different informations 
        
    This step download all the sequences in the gp format and delete it after

* Mark the position of the different proteins in the proteic sequences and the nucleic sequences.

    This step used one sequences as a representant for the other studies organisms. **It does not** seperate each sequences according to the position of their own proteins.

* Seperate the different proteins and save them 
    
    One file for each proteins.
* Identified the best model (using iqtree2)

    iqtree2 was used because it can provide a test to identify the best model for the nucleic sequences and the proteic sequences. If you only want to study the proteic sequences you may used instead prottest.

* Create the phylogenetic tree using the best model identified (using iqtree2)

* Create the caption for the Serogroup, Geographical region, and Host.

    The caption are created based on the Itol format.
    
    One color will be associated with the same values for each studies.

    The Serogroup captions are created using the **Serogroup.txt** file.

## Instructions to use 

## Instructions to adapt the code to your needs
This scripts is currently working for the peribunyaviridae (96 species) the nairoviridae(43 species) the CCHFV strains(123 strains) and the oropouche strains( strains) if you want to use this script or add other organisms you need to modifiy the following points : 

### Path and folders names
You need to change the original path and also the folders names.

To do so all you have to do is to change the variable *path* and the variables in the list *studies* in the **execution.sh** file, **download_sequences.sh**, **genbank_extract.py**, all the **mark_prot*.py** files, all the **split_protein_*.py** files and in the **All_functions.py** file.
(All the variables are at the begining of most scripts)
### E-mail adresse 

Requested are send to the NCBI databases, in order to execute these requests you need to enter your e-mail adresse.

**The package won't work if you do not enter this information**

### Add species
If you want to add another species add the name to that species in the *studies* list in the file **execution.sh** and the **download_sequences.sh**.

You have to duplicate a **mark_prot_*.py** file, **split_prot_*.py** file and **trimal_*.sh** according to your need.
(cf the description of each files to choose wich file to duplicate)

#### Mark proteins

For this script it is important to change the different position of each sequences, change the target ID according to your studies and change the startseq.

For the sutides that concern studies about different species : 
the first nine caracteres in  minusucule {organism name+proteic sequence}

The startseq was found by using blast and align the proteic sequences to the nucleic sequences. 



**If you want to add another species it is recommanded to rerun the entire script in order to have the same caption for all your studied species**

### If you you are not using the serogroup 
### If you want to use these script for other species
### If you only want to create the itol legend
### If you only want to create the Metadata files 
## Pré-requis
python 
iqtree2
package python 
phyml 
## Installation 
## Output files

### Metadata_files.txt
This file contains the following informations(if they are available in the NCBI genbank file format and the gb file format) on all the organims you choose.

        Virus name	ID itol tree	Taxonomy	Serogroup	Protein accesion number	Size Proteique sequence	Start MLD	End MLD	Sart GP38	End GP38	Start Gn	End Gn	Start NSm	End NSm	Start Gc	End Gc	Nucleique accession number	Size Nucelique sequence	Strand	Strain	Geographical region	Host	Host source	References

