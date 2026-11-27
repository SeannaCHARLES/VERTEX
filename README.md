# VERTEX

# (Viral sequence Extraction, Reconstruction and Tree EXploration)

# Project description {#sec-project-description}

This pipeline is designed to automatically download and extract protein sequences or functional domains of viral sequences from Genbank. It enables the generation of phylogenetic trees in Newick format, for each of the protein or protein domains extracted. The pipeline relies on the NCBI database and is therefore dependent on both the accuracy and the uniformity of the data provided in the input sequences. The pipeline produces iTOL-compatible captions that can be directly uploaded to the iTOL website for interactive tree visualization and additional files (e.g., FASTA alignments, metadata tables) to support downstream analyses.

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


## Launch the pipeline {#sec-launch-the-pipeline}

Make sure that the script is an executable file.

chmod -ux

When in your working directory

Enter **./\<path\>/VERTEX/Scripts/execution.sh** in your terminal.

Where \<path\> is the absolute path to the VERTEX folder that you have downloaded and decompressed from the [GitHub](https://github.com/SeannaCHARLES/VERTEX/tree/main).

## File paths {#sec-file-paths}

**these files must be in the working directory**

protein_equivalent.txt

species_equivalent.txt

all_information.txt

\<IDs_file_name\>.txt

user_legend.txt

## Pipeline description {#sec-pipeline-description}

### Execution.sh {#sec-execution.sh}

- Download all the sequences in the Genbank format.

This step performs the download of all the sequences in a Genbank format and the concatenation of all the Genbank files into one file called **sequences\<Study_name\>.gb**.

| Input | Output |
|:----------------------------|:------------------------------------------|
| all_information.txt**(From user)** | *`<Project absolute path>`*sequences\<Study_name\>.gb |

- Creates a Blast database of all the reference virus proteins.

| Input | Output | File description |
|:---------------------:|:----------------------:|-------------------------|
| \<Id_TargetID\>\_TargetID.fasta (From Genbank_extract.py) | `<Project absolute path>`protein_prot_1-n_TargetID.fasta | Fasta file with all the studied proteins of the \<TargetID**\>** |

- Launches of TrimAl

- Launches of IQtree2

Identifies the best model (using IQtree2). IQtree2 was used because it can provide a test to identify the best model of tree reconstruction for the nucleic and protein sequences. If you only want to study the protein sequences you may use prottest instead.

- Creates the phylogenetic tree using the best model identified (using IQtree2)

| Input | Output | File description |
|:-----------------:|:------------------:|----------------------------------|
| *_prot_OG.aln-fasta | `<Project absolute path>/<Res folder>/treefolder/`**<protein name>_<Study name>_prot_trimal.fasta.treefile** | Phylogenetic tree file in newick format. One file for each protein |
| protein_prot_spe.aln-fasta | *`<Project absolute path><Res folder>/treefolder/`***protein_prot_spe.treefile** | Phylogenetic tree file in newick format for the specific protein. (if you filled up the section [Specific protein](#sec-g-specific-protein)) |

- Executes of all the subsequent scripts
  
### Genbank_extract.py {#sec-genbank_extract.py}

- Separates all the nucleic sequences and the protein sequences into single sequence files. The nucleic and proteic sequences are extracted from the **sequences\<Study_name\>.gb** file and saved into 2 folders **nucl** and **prot**.

- Downloads and concatenates the GenPept(Gp) file for each sequence in the ***prot_all\_\<Study_name\>.gp*** file.

The Gp file format is very similar to the Genbank file format, except that Gp file includes the positions of the various proteins and protein domains. This format can be downloaded in the same way as the Genbank file from the protein identifier record in the NCBI database or website.

- Creates a **Metadata_files.txt** that contains different information available in the Genbank and the Gp files.

This file contains the following information:

Study \| Virus name ID for Itol tree \| Taxonomy \| Information \| Protein accession number \| Size Protein sequence \| Start proteins \| End proteins \| extract_size \| Nucleic accession number \| Size Nucleic sequence \| Strand \| Strain \| Geographical region \| Host or vector \| References

- Creates a file **pos_prot\_\<Study_name\>.txt** that lists all the positions for each protein or domain available in the Gp file.

- Creates the caption for the user_legend, Geographical region, and Host.

The captions are created based on the Itol format. One colour will be associated with the same values for all the studies.

The user-information caption is created using the **user_legend.txt** file. You can put in that file all the information you want to use to characterize your study (for example serogroup or precise taxonomy)

| Input | Output | File description |
|:------------------:|:----------------------:|----------------------------|
| sequences\<Study_name\>.gb | *`<Project absolute path>`*\<ID_nucl\>\_prot.fasta | \<ID_nucl\> protein sequence |
| user_legend.txt **(from user)** | *`<Project absolute path>`*\<ID_nucl\>\_nucl.fasta | \<ID_nucl\> nucleic sequence |
|  | `<Project absolute path>`prot_all\_.fasta | All protein sequences |
|  | `<Project absolute path>`All_nucleic_sequences.fasta | All nucleic sequences |
|  | `<Project absolute path>`All_protein_sequences_align.aln-fasta | Alignment of all protein sequences |
|  | `<Project absolute path>`All_nucleic_sequences_align.aln-fasta | Alignment of all nucleic sequences |
|  | `<Project absolute path>`All_protein_sequences.gp | Concatenated gp files. |
|  | `<Project absolute path>legend/`**Informations_caption.txt** | User caption in Itol format |
|  | `<Project absolute path>legend/`**Geographical_caption.txt** | Geographical caption in Itol format (geographical region from where the virus were isolated) |
|  | `<Project absolute path>legend/`**Host_caption.txt** | Host/vector caption in Itol format (order of the host/vector) |
|  | `<Project absolute path>legend/`**Label_caption.txt** | Virus name caption in Itol format. Changes the Itol ID into full virus name. |
|  | `<Project absolute path>legend/`**Metadata.txt** | File containing all the extracted metadata. See [genbank extract](#sec-genbank_extract.py) |
|  | `<Project absolute path>`**pos_prot\_\<Study_name\>.txt** | All the positions in the Gp files for each virus and each studied protein. |

### Prot_position.py (optional) {#sec-prot_position.py-optional}

- This file uses BlastP to complete the missing position within a study. It ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file.

| Input | Output | File description |
|:-----------------:|:------------------:|----------------------------------|
| **pos_prot\_\<Study_name\>.txt** | `<Project absolute path>`**pos_prot_final\<Study_name\>.txt** | All the positions found in the Gp files and the estimate positions using blast for each virus and each studied protein. If you have indicated "Yes" in [Blast?](#sec-o-blast) in the all_information.txt file |
| protein_prot_1-n_TargetID.fasta (from Execution.sh) |  |  |

### DNA_dist.sh {#sec-dna_dist.sh}

- Launches DNAdist from phylip and creates a divergence sequence table (`<Project absolute path>`\<Study_name\>_dnadist.txt)

| Input | Output | File description |
|:----------------------:|:--------------------:|--------------------------|
| All \_nucleic_sequences_align.aln-fasta | `<Project absolute path>`**\<Study_name\>_dnadist.txt** | A Kimura distance matrix computed from the alignment of all the nucleic sequences |

### Delete_dbl.py {#sec-delete_dbl.py}

- Deletes the sequences that are considered too closely related or too similar using the `<Project absolute path>`\<Study_name\>_dnadist.txt output table and the user defined threshold from the **all_information.txt** file.

| Input | Output | File description |
|:----------------:|:----------------:|------------------------------------|
| \<Study_name\>_dnadist.txt | `<Project absolute path>`All_nucleic_sequences_sdbl.fasta | All the nucleic sequences without the removed sequences |
|  | `<Project absolute path>`All_protein_sequences_sdbl.fasta | All the protein sequences without the removed sequences |

### Sep_prot.py {#sec-sep_prot.py}

- Separates the different proteins and saves them into individual files for each protein.

| Input | Output | File description |
|:-----------------:|:--------------------:|-------------------------------|
| All_nucleic_sequences_sdbl.fasta | `<Project absolute path><res_phylo>`*_<Study name>_nucl.fasta | All the nucleic sequences after the separation of the proteins (one file for each protein) |
| All_protein_sequences_sdbl.fasta  | `<Project absolute path><res_phylo>`*_<Study name>_prot.fasta | All the protein sequences after the separation of the proteins (one file for each protein) |
|  | `<Project absolute path><Result folder>`*_nucl_OG.fasta | All the nucleic sequences after the separation of the proteins and the chosen outgroup |
|  | `<Project absolute path><Result folder>`*_prot_OG.fasta | All the protein sequences after the separation of the proteins chosen outgroup |

### Merge_spe_prot.py (optional) {#sec-merge_spe_prot.py-optional}

- Concatenates the outgroups with their corresponding proteins, if you have chosen to use different outgroups for the proteins studied.

- Concatenates the two proteins of which you want to study the evolutionary history, if applicable.

- Aligns all the previous files.

| Input | Output |   |
|:----------------:|:----------------:|------------------------------------|
| *_<Study name>_nucl.fasta | *`<Project absolute path><Res folder>`*Nucleic_<Name of the concatenated proteins>.fasta | All the nucleic sequences of the chosen proteins to create a concatenated sequence (if you filled up the section [Specific protein](#sec-g-specific-protein)) |
| *_<Study name>_prot.fasta | *`<Project absolute path><Res folder>`*protein_<Name of the concatenated proteins>.fasta | All the protein sequences of the chosen protein to create a concatenated sequence (if you filled up the section [Specific protein](#sec-g-specific-protein)) |
| *_nucl_OG.fasta | *`<Project absolute path><Res folder>`*nucl_<Name of the concatenated proteins>.aln-fasta | Alignment of all the nucleic sequences of concatenated proteins |
| *_prot_OG.fasta | *`<Project absolute path><Res folder>`*protein_<Name of the concatenated proteins>.aln-fasta | Alignment of all the protein sequences of concatenated proteins |

### Protdist.py {#sec-protdist.py}

- Launches Protdist from phylip and create a sequence divergence table for each protein.

| Input | Output | File description |
|:-----------------:|:--------------------:|--------------------------------|
|  | *`<Project absolute path>`<Study name>_protdist_*.txt | Distance matrix table for each protein |
| protein_prot_1-n.fasta |  |  |
|  |  |  |
| protein_prot_spe.fasta |  |  |
| Label_caption.txt |  |  |

This step is executed at the end of the pipeline because it can be very long to run. It enables a pairwise comparison of each virus for each protein or protein domain, and of specific proteins whenever applicable.

## Mandatory input files {#sec-mandatory-input-files}

### all_information.txt {#sec-example-of-use}

**all_information.txt** **cannot be empty.**

This file contains all the information required for the pipeline to run correctly. It is a text file separated with **"," not tabulation**.

This file is split into two parts.

#### paths and e-mail address

The first seven lines are preceded by ‘#’; **it is important not to remove them.**

``` text
#All the information will be used in the script .......please read the read_me file for further information                                                                
#Path,<Project Absolute path>
#res_file,/<result_folder>/
#path_script,<VERTEX Absolute path>VERTEX/Script/
#type_seq,nucl;prot
#Email_adress=,<your email address>
#Name of theproject,OG_Project,ID file,Outgroup Spe (Y/N),For wich prot,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
```

The first line is simply a reminder to read the README and is not essential.

``` text
#All the information will be used in the script .......please read the read_me file for further informations
```

**The following five lines must be completed.**

``` text
#Path,<Project Absolute path>
#res_file,/<result_folder>/
#path_script,<VERTEX Absolute path>VERTEX/Script/
#type_seq,nucl;prot
#Email_adress=,<your email address>
```

The first column indicates what information is required in the second column.

The rows are organised as follows:

\<information_ask\>,\<your_value\>

You must not enter \<your_value\> in the first column. The pipeline only extracts information from the second column for all rows beginning with ‘#’.

``` text
#Path,<Project Absolute path>
```
You must enter the absolute path to the folder (working directory) in which the project folders will be created.

``` text
#res_file,/<result_folder>/
```
You must specify a folder name surrounded by ‘/’ which will be used to save the output files after the protein separation. 
No absolute path is needed

``` text
#path_script,<VERTEX Absolute path>VERTEX/Script/ 
```
You must specify the **absolute path** to the VERTEX pipeline **script folder** that you downloaded from the [GitHub](https://github.com/SeannaCHARLES/VERTEX/tree/main).

``` text
#type_seq,nucl;prot
```
You must specify the type of sequence you wish to analyse, either ‘nucl’ or ‘prot’, or both; the order does not matter.

``` text
#Email_adress=,<your email address>
```
You must enter your email address. It is essential for submitting requests to the NCBI database.

#### Your study informations

The seventh row contains headers that explain what to write in each column.

``` text
#Name of the project,Species or family,List of ID,Outgroup Spe (Y/N),for which protein,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
```

We provide here a step by step running example for the study of Retroviruses.

##### Name of the project {#sec-a-execution-project-name}

This name will be used as the working directory. It also allows you to specify which project each outgroup is associated with when several outgroups are specified (if necessary).

If you wish to analyse several segments or polyproteins from the same viral family or species, you must enter different project names; otherwise, all the files produced by the pipeline for all projects will be stored in the same folder.

**Name of the project =\> \<Project name\>**

**format : [string] no space**

``` text
Name of the project => Retrovirus_EL_Rtarticle_25
```

##### OG_Project {#sec-b-study-name}

“OG_Project” put the name of the outgroup if you want to use different outgroups for all or a few proteins.

If this does not apply to you, you can add the outgroup ID(s) to the list of studied organisms and enter the project name (cf.[Project name](#sec-a-execution-project-name)).

**OG_Project** **=\> \<Project name or OG name\>**

**format : [string] no space**

``` text
OG_Project => Retrovirus_EL_Rtarticle_25
```

##### ID file {#sec-c-id-file}

The name of the file that contains the virus IDs. **This file must be at the root of the working directory**.

**ID file** **=\> \<name of the IDs file\>**

**format : [string] no space, no absolute path**

``` text
ID file => Retrovirus_EL_Rtarticle.txt
```

##### Outgroup {#sec-d-outgroup-choice-yn}

Allows you to specify whether these viruses should be treated as outgroups or not.

**Outgroup** **=\> \<Y/N\>**

**format : [String] Yes or Y or No or N, no space**

``` text
Outgroup => N
```

##### For which protein {#sec-e-for-which-protein}

Which protein do you want to study, or for which protein will this outgroup be used.

**For which protein** **=\> \<prot1;...;protn\>**

**format : [string] separate with ";" no space**

``` text
For which protein => gag;pol;rvt;rnase;integrase
```

##### Protein names (in order) {#sec-f-protein-names-in-order}

Protein names of the [OG_Project](#sec-b-study-name)

``` text
protein name (in order) => gag;pol;rvt;rnase;integrase
```

**protein name (in order)** **=\> \<prot1;...;protn\>**

**format : [string] separate with ";" no space**

##### proteins for concatenation {#sec-g-specific-protein}

Indicates which proteins (specific proteins) are involved in case of internal viral homology, if applicable. For example, if you want to study the evolutionary history of two homologous proteins present in the same virus. This part can be empty.

**spe prot** **=\> \<speprot=prot1+prot2\>**

**format : [string] name of the specific proteins for concatenation = Protein1 +Protein2, no space**

``` text
spe prot => ''
```

##### TargetID {#sec-h-targetid}

Indicates a reference ID that will be used for BLAST alignment in case some positions are missing.

**TargetID** **=\> \<ID\>**

**format : [string] One Genbank ID, no space**

``` text
TargetID => NC_001802.1
```

##### Target/OG Proteins position {#sec-i-reference-protein-positions}

The positions of the proteins for the TargetID or for the outgroup.

**Target/OG Proteins position** **=\> \<posprot1;...;posprotn\>**

**format : [integer] positions of each protein, separate with ";" ,no space**

``` text
Target/OG Proteins position => 2;429;605;1006;1024;1144;1157;1415
```

##### Start of the proteins in the nucleic sequence {#sec-j-start-of-the-proteins-in-the-nucleic-sequence}

This position will be used as the start of the protein sequence within the nucleic sequence (for all sequences).

**Target/OG Proteins position** **=\> \<start\>**

**format : [integer] ,no space**

``` text
Start of the proteins in the nucleic sequence => 1631
```

##### DNAdist Threshold {#sec-k-dnadist-threshold}

Which threshold for the Kimura distance do you want to use.

**DNAdist Threshold =\> \<Threshold\>**

**format : [float] decimal separation with "." ,no space**

``` text
DNAdist Threshold => 0.25
```

##### Organism to keep ID {#sec-l-organism-to-keep-id}

A list of the IDs that you **do not want** to be deleted during the elimination of redundancy (too closely related sequences). It can be empty otherwise. 

**format :**

**Organism to keep ID =\> \<ID1;...;IDn\>**

**[string] separate with ";"**

``` text
Organism to keep ID => ''
```

##### Organism without a particular protein ID {#sec-m-organism-without-a-particular-protein-id}

This information is useful when certain virus do not have some proteins. You can indicate the virus ID lacking a protein and the name of the protein.

**format :**

IDs where the protein is lacking

**OUT\*\<prot\>=ID1;...;ID2**

or

IDs where the protein is present

**IN\*\<prot\>=ID1;...;ID2**

otherwise

**NONE\*NN=NN**

**[string] separate with ";" ,no space**

``` text
Organism without a particular protein ID => NONE*NN=NN
```

##### Trimal {#sec-n-trimal}

Do you want to use TrimAl for this study (Yes or No).

**format**

**Trimal =\> \<Y/N\>**

**[String] Yes or Y or No or N, no space**

``` text
Trimal => Yes
```

##### Blast {#sec-o-blast}

Do you want to use blast in case of missing positions (Yes or No).

**format**

**Blast =\> \<Y/N\>**

**[String] Yes or Y or No or N, no space**

``` text
Blast => No
```

---
In the end, your file should look like this when correctly filled:

``` text
#Name of theproject,OG_Project,ID file,Outgroup Spe (Y/N),For which protein,proteins name (in order),spe prot,TargetID,Target/OG Proteins position,Start of the proteins in the nucleic sequence,DNAdist Threshold,Organism to keep ID,ORGANISM WITH(OUT) A PARTICULAR PROTEIN ID,Trimal,blast
Retrovirus_EL_Rtarticle_25,Retrovirus_EL_Rtarticle_25,Retrovirus_EL_Rtarticle.txt,N,gag;pol;rvt;rnase;integrase,gag;pol;rvt;rnase;integrase,,NC_001802.1,2;429;605;1006;1024;1144;1157;1415,1631,0.25,,NONE*NN=NN,Yes,No
```

### Viruses IDs List (**\<IDs_file_name\>.txt)**

This file must contain the list of Genbank IDs of the studied sequences. This file is mandatory and must be in the working directory. 

``` text
AF033808.1
AF228552.1
M23385.1
M80216.1
NC_000858.1
AF033811.1
AF151794.2
GU980187.1
KT724050.1
```

### **User_legend.txt** (can be empty)

This file contains the user-provided captions. It is mandatory but will not be created by the pipeline but the file can be empty.
It is a text file separated with "\|" with a header.

``` text
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
```

### **species_equivalent.txt** (can be empty)

This file is used to modify the names of hosts/vectors that are not recognised by the NCBI databases. It is a tabular file with, in the first column, the name as it appears in the GenBank file, and in the second column, the scientific name recognised by the NCBI databases.

See [host vector section](#sec-hostvectors-caption) for further information.

``` text
#Species    Scientific Species name
Culex   Diptera
Oryzomys    Rodentia
Culicidae   Diptera
Oxymycterus Rodentia
Porcine Artiodactyla
```

### **protein_equivalent.txt**

This file is used to enable the pipeline to recognise the various names of the proteins contained in the NCBI databases

(See [section region name](#sec-region-name-in-the-database) for further information)

This is a tabular file with, in the first column, the name of the protein as it appears in the GenBank file, and in the second column, the name of that same protein as you wrote it in the section (cf.[protein names in order](#sec-f-protein-names-in-order "protein names in order")).

``` text
#Protein Database protein name
reversetranscriptase    rvt
in_ integrase
rt_rnase    rnase
gag-pro-pol gag+pol
rve integrase
rp_rtvl_h_like  rvt
rt_zfrev_like   rvt
rvt_thumb   rvt 
rt_rtv  rvt
rt_zfrev_like   rvt
```


## Skipped steps {#sec-s-steps}

If you do not want to run certain steps, you must comment out the section for these steps (indicated by comment lines in **execution.sh**). Please note that some sections must not be commented out. Also ensure that the output files for the steps you are skipping are present.

## Additional Instructions and warnings {#sec-additional-instructions-and-warnings}

### Protein positions {#sec-proteins-positions}

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

27 November 2026

# Contact {#sec-contact}
emmanuelle.lerat@univ-lyon1.fr
seanna.charles@univ-lyon1.fr
