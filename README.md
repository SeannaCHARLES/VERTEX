# VERTEX

# (Viral sequence Extraction, Reconstruction and Tree EXploration)

# Requirements {#sec-requirements}

## Bash requirements {#sec-bash-requirements}

python3 3.12.3-0ubuntu2.1

bioawk 1.0-4

iqtree2 2.0.7

seqkit 2.3.1+ds-2ubuntu0.3

trimal v1.5.rev0 build[2024-05-27]

mafft 7.505-1

phylip 1:3.697+dfsg-4build1

ncbi-blast+ 2.12.0+ds-4build2

## Python requirements {#sec-python-requirements}

Bio==1.8.4

geopy==2.5.0

matplotlib==3.11.2

numpy==2.5.3

pandas==3.0.6

Requests==2.34.2

[requirements for python](requirements.txt) & [requirements for bash](requirements_bash.txt)

![VERTEX Workflow](images/Workflow_V8.pdf)

## Mandatory input files {#sec-mandatory-input-files}

### All_information.txt {#sec-example-of-use}

**all_information.txt** **cannot be empty.**

This file contains all the information required for the pipeline to run correctly.

**all_information.txt** is a text file separated with **"," not tabulation**.

This file is split into two parts.

#### path and e-mail adresse

The first seven lines are preceded by ‘#’; **it is important not to remove them.**

---
#All the information will be used in the script .......please read the read_me file for further informations                                                                
#Path,<Project Absolute path>
#res_file,/<result_folder>/
#path_script,<VERTEX Absolute path>VERTEX/Script/
#type_seq,nucl;prot
#Email_adress=,<your email adresse>
#Name of theproject,OG_Project,ID file,Outgroup Spe (Y/N),For wich prot,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
---

``` texinfo
#All the information will be used in the script .......please read the read_me file for further informations                                                                
#Path,<Project Absolute path>
#res_file,/<result_folder>/
#path_script,<VERTEX Absolute path>VERTEX/Script/
#type_seq,nucl;prot
#Email_adress=,<your email adresse>
#Name of theproject,OG_Project,ID file,Outgroup Spe (Y/N),For wich prot,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
```

The first line is simply a reminder to read the README and is not essential.

---
#All the information will be used in the script .......please read the read_me file for further informations
---

The following five lines must be completed.

---
#Path,<Project Absolute path>
#res_file,/<result_folder>/
#path_script,<VERTEX Absolute path>VERTEX/Script/
#type_seq,nucl;prot
#Email_adress=,<your email adresse>
---

The first column indicates what information is required in the second column.

The rows are organized as follows:

\<information_ask\>,\<your_value\>

You must not enter \<your_value\> in the first column. The pipeline only extracts information from the second column for all rows beginning with ‘#’.

First line :

---
#All the information will be used in the script .......please read the read_me file for further informations
---

Second line :

---
#Path,<Project Absolute path>
---

You must enter the absolute path to the folder in which the project folders will be created.

Third line :

---
#res_file,/<result_folder>/
---

You must specify a folder name surrounded by ‘/’ which will be used to save the output files after the protein separation. (This will make it easier for you to find your files later.)

No absolute path needed

Fourth line :

---
#path_script,<VERTEX Absolute path>VERTEX/Script/ 
---

You must specify the absolute path to the VERTEX pipeline script folder that you downloaded from the [GitHub](https://github.com/SeannaCHARLES/VERTEX/tree/main).

Fifth line :

---
#type_seq,nucl;prot
---

You must specify the type of sequence you wish to analyse, either ‘nucl’ or ‘prot’, or both; the order does not matter.Sixth line :

---
#Email_adress=,<your email adresse>
---

You must enter your email address. It is essential for submitting requests to the NCBI database.

#### Your study informations

The seventh row contains headers that explain what to write in each column.

---
#Name of the project,Species or family,List of ID,Outgroup Spe (Y/N),for quel prot,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
---

We provide here a step by step running example for the study of Orthobunyaviridae with 2 differents outgroups.

##### A: Name of the project {#sec-a-execution-project-name}

This name will be used as the folder name.

It also allows you to specify which project each outgroup is associated with (if necessary).

If you wish to analyse several segments or polyproteins from the same viral family or species, you must enter two different project names; otherwise, all the files produced by the pipeline for both projects will be stored in the same folder.

---
#Name of the project,
Orthobunyaviridae,
Orthobunyaviridae,
Orthobunyaviridae,
---

##### B: OG_Project {#sec-b-study-name}

“OG_Project” This column is useful for differentiated outgroups from the list of identifiers for the viral species you wish to study.

If this does not apply to you, you can add the outgroup IDs to the list of study identifiers and enter the project name (cf.[Project name](#sec-a-execution-project-name)).

---
#Name of the project,OG_Project,
Orthobunyaviridae,Khurdun,
Orthobunyaviridae,Pacuvirus,
Orthobunyaviridae,Orthobunyaviridae,
---

##### C: ID file {#sec-c-id-file}

The name of the file that contains the virus IDs. This file must be at your root.

See

---
#Name of the project,OG_Project,ID file,
Orthobunyaviridae,Khurdun,Khurdun_file.txt,
Orthobunyaviridae,Pacuvirus,Pacuvirus_file.txt,
Orthobunyaviridae,Orthobunyaviridae,Orthobunyaviridae_file.txt,
---

##### D: Outgroup(Y/N) {#sec-d-outgroup-choice-yn}

Allows you to specify whether these viruses should be treated as outgroups or not.

---
#Name of the project,OG_Project,ID file,Outgroup,
Orthobunyaviridae,Khurdun,Khurdun_file.txt,Y,
Orthobunyaviridae,Pacuvirus,Pacuvirus_file.txt,Y,
Orthobunyaviridae,Orthobunyaviridae,Orthobunyaviridae_file.txt,N,
---

##### E: For which protein {#sec-e-for-which-protein}

Which protein do you want to study, or for which protein will this outgroup be used?

---
#Name of the project,OG_Project,ID file,Outgroup,For wich prot,
Orthobunyaviridae,Khurdun,Khurdun_file.txt,Y,nsmgn,
Orthobunyaviridae,Pacuvirus,Pacuvirus_file.txt,Y,polyprotein;gn;nsm;gc,
Orthobunyaviridae,Orthobunyaviridae,Orthobunyaviridae_file.txt,N,polyprotein;gn;nsm;gc,
---

##### F: Protein names (in order) {#sec-f-protein-names-in-order}

Protein names of the [OG_Project](#sec-b-study-name)

---
#...,ID file,Outgroup,For wich prot,
...,Khurdun_file.txt,Y,nsmgn,polyprotein;gn;gc
...,Pacuvirus_file.txt,Y,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc
...,Orthobunyaviridae_file.txt,N,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc
---

##### G: specific protein {#sec-g-specific-protein}

Indicate which proteins are involved in the study of a specific protein, if applicable.

For example, if you want to study the evolutionary history of two proteins present in the same virus.

The protein Gn of Khurdun will be used as the outgroup.

---
#...,ID file,Outgroup,For wich prot,spe prot,
...,Khurdun_file.txt,Y,nsmgn,polyprotein;gn;gc,nsmgn=gn,
...,Pacuvirus_file.txt,Y,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,,
...,Orthobunyaviridae_file.txt,N,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,nsmgn=nsm+gn
---

##### H: TargetID {#sec-h-targetid}

Enter a reference ID that will be used for BLAST alignment in case some positions are missing.

---
#...,,For wich prot,spe prot,TargetID,
...,nsmgn,polyprotein;gn;gc,nsmgn=gn,,
...,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,,,
...,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,nsmgn=nsm+gn,NC_001926.1
---

##### I: Target/OG Proteins position {#sec-i-reference-protein-positions}

The positions of the proteins for the reference or the proteins outgroup.

---
#...,For wich prot,spe prot,TargetID,Target/OG Proteins position,
...,nsmgn,polyprotein;gn;gc,nsmgn=gn,,23;292;526;944,
...,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,,,1;1458;23;305;308;476;598;1386,
...,polyprotein;gn;nsm;gc,polyprotein;gn;nsm;gc,nsmgn=nsm+gn,NC_001926.1,22;302;303;476;509;1372,
---

##### J: Start of the proteins in the nucleic sequence {#sec-j-start-of-the-proteins-in-the-nucleic-sequence}

Start of the proteins in the nucleic sequence

This position will be used as the start of the protein sequence within the nucleic acid sequence (for all sequences).

---
#...,Target/OG Proteins position,Start of the proteins in the nucleic sequence,
...,23;292;526;944,34,
...,1;1458;23;305;308;476;598;1386,53,
...,22;302;303;476;509;1372,57,
---

##### K: DNAdist Threshold {#sec-k-dnadist-threshold}

Which threshold for the kimura distance do you want to use?

---
#...,Start of the proteins in the nucleic sequence,DNAdist Threshold,
...,34,,
...,53,,
...57,0.25,
---

DNAdist Threshold

##### L: Organism to keep ID {#sec-l-organism-to-keep-id}

A list of the IDs that you **do not want** to be deleted.

**SEPARATED BY ‘;’**

---
#...,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,
...,34,,,
...,53,,,
...57,0.25,MH484340.1;NC_005775;EF990927.1
---

ID Organism to keep

##### M: Organism without a particular protein ID {#sec-m-organism-without-a-particular-protein-id}

---
#...,DNAdist Threshold,Organism to keep ID,Organism without a particular protein ID
...,,,OUT*NSM=NC_078560.1,
...,,,NONE*NN=NN,
...,0.25,MH484340.1;NC_005775;EF990927.1,NONE*NN=NN,
---

Organism without a particular protein ID

OUT\*\<prot\>=list of IDs separated by ";"

IN\*\<prot\>=liste of IDs separated by ";"

##### N: Trimal {#sec-n-trimal}

Do you want to use TrimAl for this study.

---
#...,DNAdist Threshold,Organism to keep ID,Organism without a particular protein ID,Trimal
...,,,OUT*NSM=NC_078560.1,Yes,
...,,,NONE*NN=NN,Yes,
...,0.25,MH484340.1;NC_005775;EF990927.1,NONE*NN=NN,Yes,
---

Do you want to use TrimAl

##### O: Blast {#sec-o-blast}

Do you want to use blast in case of missing positions.

---
#...,Organism without a particular protein ID,Trimal,Blast
...,OUT*NSM=NC_078560.1,Yes,Yes
...,NONE*NN=NN,Yes,Yes
...,NONE*NN=NN,Yes,Yes
---

il faut nommer clairement les ifnormations que doit contenir ce fichier en les listant ainsi que le format que ça doit avoir; le fichier exemple n'est pas simple à comprendre, un utilisateur doit pouvoir construire ce fichier là de zero; il faut aussi indiquer à quel endroit doit se trouver le fichier (et les autres) et indiquer si il faut creer un dossier qui recevra les outputs et où il doit se trouver (j'ai un probleme à se niveau pour faire tourner; j'ai l'impression que la connexion à genbank fonctionne mais où sont les fichiers récupérés ? ils ne sont pas la ou le script les attend et je ne vois pas de dossier qui les contiendrait###

### Viruses IDs List (**\<IDs_file_name\>.txt)**

This file must contain the list of Genbank IDs of the studied sequences. Ce fichier est essentiel au bon fonctionnement du pipeline.C'est une liste d'Ids. en colonne (un ids par ligne) sans entête.

Ce fichier must be at the root.

---
AF033808.1
AF228552.1
M23385.1
M80216.1
NC_000858.1
AF033811.1
AF151794.2
GU980187.1
KT724050.1
---

### **User_legend.txt**

This file contains the user-provided captions.

This file must be present and will not be created by the pipeline but it can be empty.

It is a text file seperated with "\|" avec une entête.

---
Virus name|Serogroup 
Bovine immunodeficiency|Lentivirus
Jembrana disease|Lentivirus
Bovine foamy|Spumavirus
equine foamy  Equus caballus|Spumavirus
Primate T-lymphotropic  1|Deltaretrovirus
Simian T-lymphotropic  2|Deltaretrovirus
Walleye dermal sarcoma|Epsilonretrovirus
Moloney murine leukemia|Gammaretrovirus
Reticuloendotheliosis|Gammaretrovirus
---

### **species_equivalent.txt** (can be empty)

This file is used to modify the names of hosts/vectors that are not recognised by the NCBI databases.

It is a tabular file with, in the first column, the name as it appears in the GenBank file, and in the second column, the scientific name recognised by the NCBI databases.

See [host vector section](#sec-hostvectors-caption) for further information

---
#Species	Scientific Species name
Culex	Diptera
Oryzomys	Rodentia
Culicidae	Diptera
Oxymycterus	Rodentia
Porcine	Artiodactyla
---

### **protein_equivalent.txt**

This file is used to enable the pipeline to recognise the various names of the proteins contained in the NCBI databases

(See [section region name](#sec-region-name-in-the-database) for further information)

This is a tabular file with, in the first column, the name of the protein as it appears in the GenBank file, and in the second column, the name of that same protein as you wrote it in the section (cf.[#sec-f-protein-names-in-order](#sec-f-protein-names-in-order)).

## Launch the pipeline {#sec-launch-the-pipeline}

S'assurer que le fichier est exécxutable

chmod -ux

Enter **./\<path\>/VERTEX/Scripts/execution.sh** in your terminal.

Where \<path\> is the absolute path to the VERTEX folder that you have downloaded and decompressed from the [GitHub](https://github.com/SeannaCHARLES/VERTEX/tree/main).

# Project description {#sec-project-description}

This pipeline is designed to automatically download and extract protein sequences or functional domains of viral sequences from Genbank. It enables the generation of phylogenetic trees in Newick format, for each of the protein or protein domains extracted. The pipeline relies on the NCBI database and is therefore dependent on both the accuracy and the uniformity of the data provided in the input sequences. The pipeline produces iTOL-compatible captions that can be directly uploaded to the iTOL website for interactive tree visualization and additional files (e.g., FASTA alignments, metadata tables) to support downstream analyses.

## Pipeline description {#sec-pipeline-description}

### Execution.sh {#sec-execution.sh}

- Download all the sequences in the Genbank format.

This step performs the download of all the sequences in a Genbank format and the concatenation of all the Genbank files into one file called **sequences\<Study_name\>.gb**.

| Input                              | Output                         |
|:-----------------------------------|:-------------------------------|
| all_information.txt**(From user)** | \<\>sequences\<Study_name\>.gb |

: Input-Output Execution.sh

### Genbank_extract.py {#sec-genbank_extract.py}

- Separate all the nucleic sequences and the protein sequences into single sequence files. The nucleic and proteic sequences are extracted from the **sequences\<Study_name\>.gb** file and saved in 2 folders **nucl** and **prot**.

The nucleic and protein sequences are extracted from the **sequences\<Study_name\>.gb** file and saved in 2 folders **nucl** and **prot**.

- The GenPept(Gp) file format for each sequence are downloaded and concatenated in the ***prot_all\_\<Study_name\>.gp*** file.

The Gp file format is very similar to the Genbank file format, except that Gp file includes the positions of the various proteins and protein domains. This format can be downloaded in the same way as the Genbank file from the protein identifier record in the NCBI database or website.

- Create a **Metadata_files.txt** that contains different information available in the Genbank and the Gp files.

This file contains the following information :

Study \| Virus name ID for Itol tree \| Taxonomy \| Information \| Protein accession number \| Size Protein sequence \| Start proteins \| End proteins \| extract_size \| Nucleic accession number \| Size Nucleic sequence \| Strand \| Strain \| Geographical region \| Host or vector \| References

- Create a file **pos_prot\_\<Study_name\>.txt** that lists all the positions for each protein or domain available in the Gp file.

- Create the caption for the user_legend, Geographical region, and Host.

The captions are created based on the Itol format. One colour will be associated with the same values for all the studies.

The user information caption is created using the **user_legend.txt** file. You can put in that file all the information you want to use to characterize your study (for example serogroup or precise taxonomy)

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| sequences\<Study_name\>.gb (From Genbank_extract.py) | \<ID_nucl\>\_prot.fasta |
| user_legend.txt **(from user)** | \<ID_nucl\>\_nucl.fasta |
|  | prot_all\_.fasta \| |
|  | All_nucleic_sequences.fasta |
|  | All_protein_sequences_align.aln-fasta |
|  | All_nucleic_sequences_align.aln-fasta |
|  | All_protein_sequences.gp |
|  | **Informations_caption.txt** |
|  | **Geographical_caption.txt** |
|  | **Host_caption.txt** |
|  | **Label_caption.txt** |
|  | **Metadata.txt** |
|  | **pos_prot\_\<Study_name\>.txt** |

: Input-Output Genbank.py

### Execution.sh

- Create a Blast database of all the reference virus proteins.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| \<Id_ref\>\_Ref.fasta (From Genbank_extract.py) | protein_prot_1-n_ref.fasta |

: Input-Output Execution.sh

### Prot_position.py (optional) {#sec-prot_position.py-optional}

- This file use blastp to complete the missing position within a study.

This file ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| **pos_prot\_\<Study_name\>.txt** (From Genbank.py) | **pos_prot_final\<Study_name\>.txt** |
| protein_prot_1-n_ref.fasta (from Execution.sh) |  |

: Input-Output Prot_position.py

### DNA_dist.sh {#sec-dna_dist.sh}

- Launch DNAdist from phylip and create a divergence sequence table.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| All \_nucleic_sequences_align.aln-fasta (From Genbank_extract.py) | **Divergence_nucleic_table.txt** |

: Input-Output DNA_dist.sh

### Delete_dbl.py {#sec-delete_dbl.py}

- Delete the sequences that are considered too closely related or too similar using the **DNA_dist.sh** output table and the user defined threshold from the **All_information.txt** file.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| Divergence_nucleic_table.txt(From DNA_dist.sh) | All_nucleic_sequences_sdbl.fasta |
|  | All_protein_sequences_sdbl.fasta |

: Input-Output Delete_dbl.py

### Sep_prot.py {#sec-sep_prot.py}

- Separate the different proteins and save them in separated files for each protein.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| All_nucleic_sequences_sdbl.fasta (From Delete_dbl.py) | Nucleic_prot_1-n.fasta |
| All_protein_sequences_sdbl.fasta (From Delete_dbl.py) | protein_prot_1-n.fasta |
|  | Nucleic_prot_OG.fasta |
|  | protein_prot_OG.fasta |

: Input-Output Sep_prot.py

### Merge_spe_prot.py (optional) {#sec-merge_spe_prot.py-optional}

- Concatenates the outgroups with their corresponding proteins, if you have chosen to use different outgroups for the proteins studied.

- Concatenates the two proteins of which you want to study the evolutionary history, if applicable.

- Aligns all the previous files.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| Nucleic_prot_1-n.fasta (From Test_sript_Sep_prot.py) | Nucleic_prot_spe.fasta |
| Protein_prot_1-n.fasta (From Test_sript_Sep_prot.py) | protein_prot_spe.fasta |
| Nucleic_prot_OG.fasta (From Test_sript_Sep_prot.py) | protein_prot_spe.aln-fasta |
| Protein_prot_OG.fasta (From Test_sript_Sep_prot.py) | protein_prot_1-n.aln-fasta |

: Input-Output Merge_spe_prot.py

### Execution.sh {#sec-execution.sh2}

- TrimAl

- IQtree2

Identify the best model (using iqtree2)

IQtree2 was used because it can provide a test to identify the best model of tree reconstruction for the nucleic and protein sequences. If you only want to study the protein sequences you may use prottest instead.

- Create the phylogenetic tree using the best model identified (using IQtree2)

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| protein_prot_1-n.aln-fasta (From Merge_spe_prot.py) | **protein_prot_1-n.treefile** |
| protein_prot_spe.aln-fasta (From Merge_spe_prot.py) | **protein_prot_spe.treefile** |

: Input-Output Execution.sh

### Protdist.py {#sec-protdist.py}

- Launch Protdist from phylip and create a sequence divergence table for each protein.

| Input | Output |
|:----------------------------------:|:----------------------------------:|
| Nucleic_prot_1-n.fasta (From Sep_prot.py) | **Divergence_protein_table.py** |
| protein_prot_1-n.fasta (From Sep_prot.py) |  |
| Nucleic_prot_spe.fasta (From Merge_spe_prot.py) |  |
| protein_prot_spe.fasta (From Merge_spe_prot.py) |  |
| Label_caption.txt (From Genbank_extract.py) |  |

: Input-Output Protdist.py

This step is executed at the end of the pipeline because it can be very long to run. It enables a pairwise comparison of each virus for each protein or protein domain, and of specific proteins whenever applicable.

# Instructions to use {#sec-instructions-to-use}

## Informations_file.txt {#sec-informations_file.txt}

This file is the execution input file. You need to fill up all the information in this file in order to use the pipeline.

**!If an information is missing it may not work!**

## Separator {#sec-separator}

DO NOT USE "," AS A SEPARATOR INSIDE EACH COLUMN USE ";" INSTEAD

## File paths {#sec-file-paths}

protein_equivalent.txt **must be at the root**

species_equivalent.txt **must be at the root**

all_information.txt **must be at the root**

\<IDs_file_name\>.txt **must be at the root**

user_legend.txt **must be at the root**

## Skip steps {#sec-s-steps}

If you do not want to run certain steps, you must comment out the section for these steps (indicated by comment lines in **Execution.sh**). Please note that some sections must not be commented out. Also ensure that the output files for the steps you are skipping are present.

## Additional Instructions and warnings {#sec-additional-instructions-and-warnings}

### Proteins positions {#sec-proteins-positions}

Please pay attention to the ***Final_pos_prot.txt*** file. The blast alignment is only an approximation Please pay attention to the ***Final_pos_prot.txt*** file. The blast alignment is only an approximation

This is not a completely reliable solution. We strongly advise against using VERTEX if you have very few virus with records in the databases.

It is also strongly recommended that you check the consistency between the actual positions of the proteins and the predicted positions. The **Pos_prot.txt** file contains all the positions extracted from the database.

If you are not satisfied with the positions, you can change them manually. After that, simply restart the pipeline from the **Sep_prot.py** file.

### Region name in the database

The same proteins may have different region names from one virus to another. If a protein returns very few results and you know this is unexpected, its name may not have been recognised. You can use the file ***prot_all\_\<Study\>.gp*** to help you identified the different region names.

Whether it differs from the name you entered in the various protein names (red rectangle).

![Gc and Gn different name](images/Diff_G1_G2.pdf)

Or that it differs from the various names of the same protein (blue rectangle)

![Polyprotein different name](images/diff_nom_polyprot.pdf)

You can add protein name equivalents to the file **protein_equivalent.txt**

such as: ***glycoproteinprecursor\<\\t\>polyprotein***

### Host/vectors caption {#sec-hostvectors-caption}

The names of the hosts/vectors are not all uniform and some are not recognised by the NCBI database. To get pass this, you can add equivalents for the host/vector name and the scientific name in the **species_equivalent.txt**, such as: ***mosquito\<\\t\> Culicidae***

The images of the various orders of hosts/vectors are available on another public GitHub repository [Host/vector images](https://github.com/SeannaCHARLES/Images "Host/vector images").

# Citations {#sec-citations}

If using this tool, please cite Seanna Charles, Anja Böckmann and Emmanuelle Lerat (2026) VERTEX: A pipeline for Viral sequence Extraction, Reconstruction and Tree EXploration.

12 November 2026

# Contact {#sec-contact}
