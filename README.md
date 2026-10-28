---
---
---

# VERTEX (Viral sequence Extraction, Reconstruction and Tree EXploration)

## requirements

### Bash requirements

python3 3.12.3-0ubuntu2.1

bioawk 1.0-4

iqtree2 2.0.7

seqkit 2.3.1+ds-2ubuntu0.3

trimal v1.5.rev0 build[2024-05-27]

mafft 7.505-1

phylip 1:3.697+dfsg-4build1

ncbi-blast+ 2.12.0+ds-4build2

### Python requirements

Bio==1.8.4

geopy==2.5.0

matplotlib==3.11.2

numpy==2.5.3

pandas==3.0.6

Requests==2.34.2

[requirements for python](requirements.txt) & [requirements for bash](requirements_bash.txt)

![VERTEX Workflow](images/Workflow_V8.pdf)

## Project description

This pipeline is designed to automatically download and extract protein sequences or functional domains of viral sequences from Genbank. It enables the generation of phylogenetic trees in Newick format, for each of the protein or protein domains extracted. The pipeline relies on the NCBI database and is therefore dependent on both the accuracy and the uniformity of the data provided in the input sequences. The pipeline produces iTOL-compatible captions that can be directly uploaded to the iTOL website for interactive tree visualization and additional files (e.g., FASTA alignments, metadata tables) to support downstream analyses.

## Pipeline description

### Execution.sh

- Download all the sequences in the Genbank format.

This step performs the download of all the sequences in a Genbank format and the concatenation of all the Genbank files into one file called **sequences\<Study_name\>.gb**.

+---------------------+----------------------------+
| Input               | Output                     |
+:===================:+:==========================:+
| all_information.txt | sequences\<Study_name\>.gb |
|                     |                            |
| **(From user)**     |                            |
+---------------------+----------------------------+

: Input-Output Execution.sh

### Genbank_extract.py

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

+----------------------------+---------------------------------------+
| Input                      | Output                                |
+:==========================:+:=====================================:+
| sequences\<Study_name\>.gb | \<ID_nucl\>\_prot.fasta               |
|                            |                                       |
| (From Genbank_extract.py)  | \<ID_nucl\>\_nucl.fasta               |
|                            |                                       |
| user_legend.txt            | prot_all\_.fasta \|                   |
|                            |                                       |
| **(from user)**            | All_nucleic_sequences.fasta           |
|                            |                                       |
|                            | All_protein_sequences_align.aln-fasta |
|                            |                                       |
|                            | All_nucleic_sequences_align.aln-fasta |
|                            |                                       |
|                            | All_protein_sequences.gp              |
|                            |                                       |
|                            | **Informations_caption.txt**          |
|                            |                                       |
|                            | **Geographical_caption.txt**          |
|                            |                                       |
|                            | **Host_caption.txt**                  |
|                            |                                       |
|                            | **Label_caption.txt**                 |
|                            |                                       |
|                            | **Metadata.txt**                      |
|                            |                                       |
|                            | **pos_prot\_\<Study_name\>.txt**      |
+----------------------------+---------------------------------------+

: Input-Output Genbank.py

### Execution.sh

- Create a Blast database of all the reference virus proteins.

+---------------------------+----------------------------+
| Input                     | Output                     |
+:=========================:+:==========================:+
| \<Id_ref\>\_Ref.fasta     | protein_prot_1-n_ref.fasta |
|                           |                            |
| (From Genbank_extract.py) |                            |
+---------------------------+----------------------------+

: Input-Output Execution.sh

### Prot_position.py (optional)

- This file use blastp to complete the missing position within a study.

This file ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file.

+----------------------------------+--------------------------------------+
| Input                            | Output                               |
+:================================:+:====================================:+
| **pos_prot\_\<Study_name\>.txt** | **pos_prot_final\<Study_name\>.txt** |
|                                  |                                      |
| (From Genbank.py)                |                                      |
|                                  |                                      |
| protein_prot_1-n_ref.fasta       |                                      |
|                                  |                                      |
| (from Execution.sh)              |                                      |
+----------------------------------+--------------------------------------+

: Input-Output Prot_position.py

### DNA_dist.sh

- Launch DNAdist from phylip and create a divergence sequence table.

+-----------------------------------------+----------------------------------+
| Input                                   | Output                           |
+:=======================================:+:================================:+
| All \_nucleic_sequences_align.aln-fasta | **Divergence_nucleic_table.txt** |
|                                         |                                  |
| (From Genbank_extract.py)               |                                  |
+-----------------------------------------+----------------------------------+

: Input-Output DNA_dist.sh

### Delete_dbl.py

- Delete the sequences that are considered too closely related or too similar using the **DNA_dist.sh** output table and the user defined threshold from the **All_information.txt** file.

+------------------------------+----------------------------------+
| Input                        | Output                           |
+:============================:+:================================:+
| Divergence_nucleic_table.txt | All_nucleic_sequences_sdbl.fasta |
|                              |                                  |
| (From DNA_dist.sh)           | All_protein_sequences_sdbl.fasta |
+------------------------------+----------------------------------+

: Input-Output Delete_dbl.py

### Sep_prot.py

- Separate the different proteins and save them in separated files for each protein.

+----------------------------------+------------------------+
| Input                            | Output                 |
+:================================:+:======================:+
| All_nucleic_sequences_sdbl.fasta | Nucleic_prot_1-n.fasta |
|                                  |                        |
| (From Delete_dbl.py)             | protein_prot_1-n.fasta |
|                                  |                        |
| All_protein_sequences_sdbl.fasta | Nucleic_prot_OG.fasta  |
|                                  |                        |
| (From Delete_dbl.py)             | protein_prot_OG.fasta  |
+----------------------------------+------------------------+

: Input-Output Sep_prot.py

### Merge_spe_prot.py (optional)

- Align all proteins with their outgroup and create a file containing the %%%proteins you wish to use for the specific protein%%% ce n'est pas tres clair

+-------------------------------+----------------------------+
| Input                         | Output                     |
+:=============================:+:==========================:+
| Nucleic_prot_1-n.fasta        | Nucleic_prot_spe.fasta     |
|                               |                            |
| (From Test_sript_Sep_prot.py) | protein_prot_spe.fasta     |
|                               |                            |
| protein_prot_1-n.fasta        | protein_prot_spe.aln-fasta |
|                               |                            |
| (From Test_sript_Sep_prot.py) | protein_prot_1-n.aln-fasta |
|                               |                            |
| Nucleic_prot_OG.fasta         |                            |
|                               |                            |
| (From Test_sript_Sep_prot.py) |                            |
|                               |                            |
| protein_prot_OG.fasta         |                            |
|                               |                            |
| (From Test_sript_Sep_prot.py) |                            |
+-------------------------------+----------------------------+

: Input-Output Merge_spe_prot.py

### Execution.sh

- TrimAl

- IQtree2

Identify the best model (using iqtree2)

IQtree2 was used because it can provide a test to identify the best model of tree reconstruction for the nucleic and protein sequences. If you only want to study the protein sequences you may use prottest instead.

- Create the phylogenetic tree using the best model identified (using IQtree2)

+----------------------------+-------------------------------+
| Input                      | Output                        |
+:==========================:+:=============================:+
| protein_prot_1-n.aln-fasta | **protein_prot_1-n.treefile** |
|                            |                               |
| (From Merge_spe_prot.py)   | **protein_prot_spe.treefile** |
|                            |                               |
| protein_prot_spe.aln-fasta |                               |
|                            |                               |
| (From Merge_spe_prot.py)   |                               |
+----------------------------+-------------------------------+

: Input-Output Execution.sh

### Protdist.py

- Launch Protdist from phylip and create a sequence divergence table for each protein.

+---------------------------+---------------------------------+
| Input                     | Output                          |
+:=========================:+:===============================:+
| Nucleic_prot_1-n.fasta    | **Divergence_protein_table.py** |
|                           |                                 |
| (From Sep_prot.py)        |                                 |
|                           |                                 |
| protein_prot_1-n.fasta    |                                 |
|                           |                                 |
| (From Sep_prot.py)        |                                 |
|                           |                                 |
| Nucleic_prot_spe.fasta    |                                 |
|                           |                                 |
| (From Merge_spe_prot.py)  |                                 |
|                           |                                 |
| protein_prot_spe.fasta    |                                 |
|                           |                                 |
| (From Merge_spe_prot.py)  |                                 |
|                           |                                 |
| Label_caption.txt         |                                 |
|                           |                                 |
| (From Genbank_extract.py) |                                 |
+---------------------------+---------------------------------+

: Input-Output Protdist.py

This step is executed at the end of the pipeline because it can be very long to run. It enables a pairwise comparison of each virus for each protein or protein domain, and of specific proteins whenever applicable.

## Instructions to use

### Informations_file.txt

This file is the execution input file. You need to fill up all the information in this file in order to use the pipeline.

**!If an information is missing it may not work!**

### Separator

DO NOT USE "," AS A SEPARATOR INSIDE EACH COLUMN USE ";" INSTEAD

### File paths

**protein_equivalent.txt** must be at the root

**species_equivalent.txt** must be at the root

**all_information.txt** must be at the root

**Virus_ids.txt** must be at the root

**user_informations.txt** must be at the root

### Example of use

We provide here a step by step running example for the study of Orthobunyaviridae with a specific outgroup

#### A: Execution Project name

This name will be used as the folder name. It also allows you to specify which project the outgroup-specific data is associated with, where necessary.

If you wish to analyse several segments or polyproteins from the same viral family or species, you must enter two different project names; otherwise, all the files produced by the pipeline for both projects will be stored in the same folder.

+-------------------+
| Project name      |
+===================+
| Orthobunyaviridae |
+-------------------+
| Orthobunyaviridae |
+-------------------+
| Orthobunyaviridae |
+-------------------+

: Project name

#### B: Study name %%cette partie n'est pas tres claire pour moi%%%

Study name permet de différencier les informations qui concernent des outgroup qui ne serait pas utilisé pour toutes les protéines de l'étude.

Par exemple si vous avez une protéine spécifique, ou si vous voulez utiliser des outgroups différents pour chaque protéine. Si ce n'est pas votre cas vous pouvez ajouter l'outgroup dans la liste des identifiants de l'étude et saisir le nom du projet.

This name is used to distinguish specific outgroups (OG).

If the row ??%%terme pas clair%% does not concern a specific outgroup, enter the project name. Otherwise, enter a name for the specific outgroup (ex. name).

If you wish to use different specific outgroups for the same project, you must have one row per outgroup.

+-------------------+-------------------+
| A                 | Study name        |
+===================+===================+
| Orthobunyaviridae | Khurdun           |
+-------------------+-------------------+
| Orthobunyaviridae | Pacuvirus         |
+-------------------+-------------------+
| Orthobunyaviridae | Orthobunyaviridae |
+-------------------+-------------------+

: Study name : Khurdun and Pacuvirus will be used as outgroups

#### C: ID file

The name of the file that contains the virus IDs. This file must be at the root.

+-------------------+-------------------+----------------------------+
| A                 | B                 | ID files                   |
+===================+===================+============================+
| Orthobunyaviridae | Khurdun           | Khurdun_file.txt           |
+-------------------+-------------------+----------------------------+
| Orthobunyaviridae | Pacuvirus         | Pacuvirus_file.txt         |
+-------------------+-------------------+----------------------------+
| Orthobunyaviridae | Orthobunyaviridae | Orthobunyaviridae_file.txt |
+-------------------+-------------------+----------------------------+

: ID file

#### D: Outgroup choice (Y/N)

Permet de spécifier si ces virus doivent être considéré comme des outgroups ou non.

%%line =\> terme pas clair%% Y : if this line concern a specific outgroup for one or multiple proteins.

N : Is it a outgroup.

+-------------------+-------------------+-----+--------------+
| A                 | B                 | C   | Outgroup Spe |
+===================+===================+=====+==============+
| Orthobunyaviridae | Khurdun           | IDk | Y            |
+-------------------+-------------------+-----+--------------+
| Orthobunyaviridae | Pacuvirus         | IDp | Y            |
+-------------------+-------------------+-----+--------------+
| Orthobunyaviridae | Orthobunyaviridae | IDn | N            |
+-------------------+-------------------+-----+--------------+

: Outgroup Specific : Is it a specifc outgroup?

#### E: For which protein

Quelle protéine voulez étudier , ou pour quelle protéine cet outgroup sera t-il utilisé.

Which proteins are concerned by this study/line?

%%pas clair pour moi, surtout ce qui est écrit en dessous%%% If this line doesn't concern a specific outgroup, enter the names of the proteins from the studied species .

If this line refers to a specific outgroups, specify for which proteins this outgroup will be used

+-------------------+-------------------+---+-----------------------+
| A                 | B                 | D | For wich prot         |
+===================+===================+===+=======================+
| Orthobunyaviridae | Khurdun           | Y | nsmgn                 |
+-------------------+-------------------+---+-----------------------+
| Orthobunyaviridae | Pacuvirus         | Y | polyprotein;gn;nsm;gc |
+-------------------+-------------------+---+-----------------------+
| Orthobunyaviridae | Orthobunyaviridae | N | polyprotein;gn;nsm;gc |
+-------------------+-------------------+---+-----------------------+

: For which proteins

#### F: Protein names (in order)

Protein names of the study name (col B)

+-------------------+-------------------+-------+--------------------------+
| A                 | B                 | E     | Proteins name (in order) |
+===================+===================+=======+==========================+
| Orthobunyaviridae | Khurdun           | nsmgn | polyprotein;gn;gc        |
+-------------------+-------------------+-------+--------------------------+
| Orthobunyaviridae | Pacuvirus         | prots | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+-------+--------------------------+
| Orthobunyaviridae | Orthobunyaviridae | prots | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+-------+--------------------------+

: Protein names (in order)

#### G: specific protein

Identify the proteins that will be involved in the specific protein %%phrase pas claire%%

+-------------------+-------------------+-------+--------------+
| A                 | B                 | F     | spe prot     |
+===================+===================+=======+==============+
| Orthobunyaviridae | Khurdun           | prots | nsmgn=gn     |
+-------------------+-------------------+-------+--------------+
| Orthobunyaviridae | Pacuvirus         | prots |              |
+-------------------+-------------------+-------+--------------+
| Orthobunyaviridae | Orthobunyaviridae | prots | nsmgn=nsm+gn |
+-------------------+-------------------+-------+--------------+

: specific protein

The protein Gn of Khurdun will be used as the outgroup.

#### H: TargetID

Enter a reference ID that will be used for BLAST alignment in case some positions are missing.

+-------------------+-------------------+--------------+-------------+
| A                 | B                 | G            | RefID       |
+===================+===================+==============+=============+
| Orthobunyaviridae | Khurdun           | nsmgn=gn     |             |
+-------------------+-------------------+--------------+-------------+
| Orthobunyaviridae | Pacuvirus         |              |             |
+-------------------+-------------------+--------------+-------------+
| Orthobunyaviridae | Orthobunyaviridae | nsmgn=nsm+gn | NC_001926.1 |
+-------------------+-------------------+--------------+-------------+

: RefID

#### I: reference protein positions

The positions of the proteins for the reference or the specific outgroup.

+-------------------+-------------------+-------+---------------------------------+
| A                 | B                 | H     | reference /OG proteins position |
+===================+===================+=======+=================================+
| Orthobunyaviridae | Khurdun           |       | 23;292;526;944                  |
+-------------------+-------------------+-------+---------------------------------+
| Orthobunyaviridae | Pacuvirus         |       | 1;1458;23;305;308;476;598;1386  |
+-------------------+-------------------+-------+---------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDref | 22;302;303;476;509;1372         |
+-------------------+-------------------+-------+---------------------------------+

: RefID

#### J: Start of the proteins in the nucleic sequence

+-------------------+-------------------+-----+-----------------------------------------------+
| A                 | B                 | I   | Start of the proteins in the nucleic sequence |
+===================+===================+=====+===============================================+
| Orthobunyaviridae | Khurdun           | pos | 34                                            |
+-------------------+-------------------+-----+-----------------------------------------------+
| Orthobunyaviridae | Pacuvirus         | pos | 53                                            |
+-------------------+-------------------+-----+-----------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | pos | 57                                            |
+-------------------+-------------------+-----+-----------------------------------------------+

: Start of the proteins in the nucleic sequence

This position will be used as the start of the protein sequence within the nucleic acid sequence (for all sequences).

**K: DNAdist Threshold**

+-------------------+-------------------+----+-------------------+
| A                 | B                 | J  | DNAdist Threshold |
+===================+===================+====+===================+
| Orthobunyaviridae | Khurdun           | 34 | 0.01              |
+-------------------+-------------------+----+-------------------+
| Orthobunyaviridae | Pacuvirus         | 53 | 0.05              |
+-------------------+-------------------+----+-------------------+
| Orthobunyaviridae | Orthobunyaviridae | 57 | 0.01              |
+-------------------+-------------------+----+-------------------+

: DNAdist Threshold

#### L: Organism to keep ID

a list of the IDs that you do not want to be deleted.

**SEPARATED BY ‘;’**

+-------------------+-------------------+------+---------------------------------+
| A                 | B                 | k    | Organism to keep ID             |
+===================+===================+======+=================================+
| Orthobunyaviridae | Khurdun           | 0.01 |                                 |
+-------------------+-------------------+------+---------------------------------+
| Orthobunyaviridae | Pacuvirus         | 0.05 |                                 |
+-------------------+-------------------+------+---------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | 0.01 | MH484340.1;NC_005775;EF990927.1 |
+-------------------+-------------------+------+---------------------------------+

: ID Organism to keep

#### M: Organism without a particular protein ID

+-------------------+-------------------+-----+------------------------------------------+
| A                 | B                 | L   | Organism without a particular protein ID |
+===================+===================+=====+==========================================+
| Orthobunyaviridae | Khurdun           |     | OUT\*NSM=NC_078560.1                     |
+-------------------+-------------------+-----+------------------------------------------+
| Orthobunyaviridae | Pacuvirus         |     | NONE\*NN=NN                              |
+-------------------+-------------------+-----+------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDs | NONE\*NN=NN                              |
+-------------------+-------------------+-----+------------------------------------------+

: Organism without a particular protein ID

OUT\*\<prot\>=list of IDs

IN\*\<prot\>=liste of IDs

**N: Trimal**

Do you want to use TrimAl for this study.

+-------------------+-------------------+-------------+--------+
| A                 | B                 | M           | Trimal |
+===================+===================+=============+========+
| Orthobunyaviridae | Khurdun           | OUT\*NSM=ID | Yes    |
+-------------------+-------------------+-------------+--------+
| Orthobunyaviridae | Pacuvirus         | NONE        | Yes    |
+-------------------+-------------------+-------------+--------+
| Orthobunyaviridae | Orthobunyaviridae | NONE        | Yes    |
+-------------------+-------------------+-------------+--------+

: Do you want to use TrimAl

**O: Blast**

Do you want to use blast in case of missing proteins positions

+-------------------+-------------------+---+-------+
| A                 | B                 | N | Blast |
+===================+===================+===+=======+
| Orthobunyaviridae | Khurdun           | Y | Y     |
+-------------------+-------------------+---+-------+
| Orthobunyaviridae | Pacuvirus         | Y | Y     |
+-------------------+-------------------+---+-------+
| Orthobunyaviridae | Orthobunyaviridae | Y | Y     |
+-------------------+-------------------+---+-------+

: Do you want to use blast?

## Additional Instructions and warnings

### Proteins positions

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

### Host/vectors caption

The names of the hosts/vectors are not all uniform and some are not recognized by the NCBI database. To get pass this, you can add equivalents for the host/vector name and the scientific name in the **species_equivalent.txt**

such as: ***mosquito\<\\t\> Culicidae***

The Images of the various orders of hosts/vectors are available on another public GitHub repository. Feel free to add new illustrations, this GitHub repository does not yet contain all the orders.

The images %%quelles images ? celles que tu mentionnes au dessus?%%% were generated by AI.

## Citation
If using this tool, please cite Seanna Charles, Anja Böckmann and Emmanuelle Lerat (2026) VERTEX: A pipeline for Viral sequence Extraction, Reconstruction and Tree EXploration.

28 October 2026
