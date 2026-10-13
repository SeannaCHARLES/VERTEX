---
---
---

# VERTEX

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

### Pyhton requirements

Bio==1.8.4

geopy==2.5.0

matplotlib==3.11.2

numpy==2.5.3

pandas==3.0.6

Requests==2.34.2

[requirements for python](requirements.txt) & [requirements for bash](requirements_bash.txt)

![Workflow of the vertex](images/Workflow_v6.pdf)

## Project description

The purpose of this pipeline is to automatically separate proteins or proteic domains from sequences, when you don't have a corresponding ID in the NCBI database. It can generate phylogenetic trees in newick format with the separate regions.

This pipeline use the NCBI database and is therefore dependent not only on the information provided but also on the uniformity of that information.

It will create caption that you can directly use on the website Itol. It will also allow you to create other files that may be useful for your research.

## Pipeline description  

### Execution.sh

- Download all the sequences in the Genbank format.

This step download all the sequences in a Genbank format and concatenated all the Genbank files into one file name **sequences\<Study_name\>.gb**.

+---------------------+----------------------------+
| Input               | Output                     |
+:===================:+:==========================:+
| all_information.txt | sequences\<Study_name\>.gb |
|                     |                            |
| **(From user)**     |                            |
+---------------------+----------------------------+

: Input-Output Execution.sh

### Genbank_extract.py

- Separate all the nucleic sequences and the proteic sequences into single sequence file.

From the **sequences\<Study_name\>.gb** file the nucleic and proteic sequences are being extract and save in 2 folders **nucl** and **prot**.

- The gp file format for each sequence are being downloaded and concatenate in the ***prot_all\_\<Study_name\>.gp*** file.
- Create a **Metadata_files.txt** that contains different information available in the genbank file format and the Gp file format.

This file contains the following information : Study–Virus name–ID itol tree–Taxonomy–Informations–Protein accession number–Size Proteic sequence–Startproteins–Endproteines–Nucleic accession number–Size Nucleic sequence–Strand–Strain–Geographical region–Host or vector References

- Create a file that lits all the positions for each protein or domain available in the gp files format saved in the file **pos_prot\_\<Study_name\>.txt**.
- Create the caption for the user_informations, Geographical region, and Host.

The caption are created based on the Itol format. One color will be associated with the same values for all the studies.

The personal information caption is created using the **user_informations.txt** file. You can put in that file all the information you want to use to caracterise you study (for example serogroup or precise taxonomy)

+----------------------------+---------------------------------------+
| Input                      | Output                                |
+:==========================:+:=====================================:+
| sequences\<Study_name\>.gb | \<ID_nucl\>\_prot.fasta               |
|                            |                                       |
| (From Genbank_extract.py)  | \<ID_nucl\>\_nucl.fasta               |
|                            |                                       |
| user_informations.txt      | prot_all\_.fasta \|                   |
|                            |                                       |
| **(from user)**            | All_nucleic_sequences.fasta           |
|                            |                                       |
|                            | All_proteic_sequences_align.aln-fasta |
|                            |                                       |
|                            | All_nucleic_sequences_align.aln-fasta |
|                            |                                       |
|                            | All_proteic_sequences.gp              |
|                            |                                       |
|                            | **Informations_caption.txt**          |
|                            |                                       |
|                            | **Geographical_caption.txt**          |
|                            |                                       |
|                            | **Host_caption.txt**                  |
|                            |                                       |
|                            | **Metadata.txt**                      |
|                            |                                       |
|                            | **pos_prot\_\<Study_name\>.txt**      |
+----------------------------+---------------------------------------+

: Input-Output Genbank.py

### Execution.sh

- Create a Blast database of all the reference virus/organisms proteins.

+----------------------------------+-----------------------------------+
| Input                            | Output                            |
+:================================:+:=================================:+
| \<Id_ref\>\_Ref.fasta            | Proteic_prot_1-n_ref.fasta        |
|                                  |                                   |
| (From Genbank_extract.py)        |                                   |
+----------------------------------+-----------------------------------+

: Input-Output Execution.sh

### Prot_position.py (optional)

- This file use blastp to complete the missing position within a study.

This file ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file format.

+----------------------------------+--------------------------------------+
| Input                            | Output                               |
+:================================:+:====================================:+
| **pos_prot\_\<Study_name\>.txt** | **pos_prot_final\<Study_name\>.txt** |
|                                  |                                      |
| (From Genbank.py)                |                                      |
|                                  |                                      |
| Proteic_prot_1-n_ref.fasta       |                                      |
|                                  |                                      |
| (from Execution.sh)              |                                      |
+----------------------------------+--------------------------------------+

: Input-Output Prot_position.py

### DNA_dist.sh

- This file lauch DNA dist from phylip and create a divergence sequence table.

+-----------------------------------------+----------------------------------+
| Input                                   | Output                           |
+:=======================================:+:================================:+
| All \_nucleic_sequences_align.aln-fasta | **Divergence_nucleic_table.txt** |
|                                         |                                  |
| (From Genbank_extract.py)               |                                  |
+-----------------------------------------+----------------------------------+

: Input-Output DNA_dist.sh

### Delete_dbl.py

- This file delete de sequences that are judged too close using the **DNA_dist.sh** output table and the define threshold define in the **All_information.txt (from user)**.

+------------------------------+----------------------------------+
| Input                        | Output                           |
+:============================:+:================================:+
| Divergence_nucleic_table.txt | All_nucleic_sequences_sdbl.fasta |
|                              |                                  |
| (From DNA_dist.sh)           | All_proteic_sequences_sdbl.fasta |
+------------------------------+----------------------------------+

: Input-Output Delete_dbl.py

### Mark_spe.py

- Separate the different proteins and save them in one file for each proteins.

+----------------------------------+------------------------+
| Input                            | Output                 |
+:================================:+:======================:+
| All_nucleic_sequences_sdbl.fasta | Nucleic_prot_1-n.fasta |
|                                  |                        |
| (From Delete_dbl.py)             | Proteic_prot_1-n.fasta |
|                                  |                        |
| All_proteic_sequences_sdbl.fasta | Nucleic_prot_OG.fasta  |
|                                  |                        |
| (From Delete_dbl.py)             | Proteic_prot_OG.fasta  |
+----------------------------------+------------------------+

: Input-Output Test_sript_mark_spe.py

### Merge_spe_prot.py (optional)

- This file merge the specific OG species, with the specific study.

+-------------------------------+------------------------+
| Input                         | Output                 |
+:=============================:+:======================:+
| Nucleic_prot_1-n.fasta        | Nucleic_prot_spe.fasta |
|                               |                        |
| (From Test_sript_mark_spe.py) | Proteic_prot_spe.fasta |
|                               |                        |
| Proteic_prot_1-n.fasta        |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
|                               |                        |
| Nucleic_prot_OG.fasta         |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
|                               |                        |
| Proteic_prot_OG.fasta         |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
+-------------------------------+------------------------+

: Input-Output Merge_spe_prot.py

### Protdist.py

- This file lauch Protdist from phyllip and create a sequence divergence table for each protein.

+------------------------------+---------------------------------+
| Input                        | Output                          |
+:============================:+:===============================:+
| Nucleic_prot_1-n.fasta       | **Divergence_proteic_table.py** |
|                              |                                 |
| (From mark_spe.py)           |                                 |
|                              |                                 |
| Proteic_prot_1-n.fasta       |                                 |
|                              |                                 |
| (From mark_spe.py)           |                                 |
|                              |                                 |
| Nucleic_prot_spe.fasta       |                                 |
|                              |                                 |
| (From Merge_spe_prot.py)     |                                 |
|                              |                                 |
| Proteic_prot_spe.fasta       |                                 |
|                              |                                 |
| (From Merge_spe_prot.py)     |                                 |
+------------------------------+---------------------------------+

: Input-Output Protdist.py

### Execution.sh

- Trimal

- IQtree2

Identified the best model (using iqtree2)

iqtree2 was used because it can provide a test to identify the best model for the nucleic sequences(not present in this pipeline) and the proteic sequences. If you only want to study the proteic sequences you may used prottest instead.

- Create the phylogenetic tree using the best model identified (using iqtree2)

+------------------------------+-------------------------------+
| Input                        | Output                        |
+:============================:+:=============================:+
| Proteic_prot_1-n.fasta       | **Proteic_prot_1-n.treefile** |
|                              |                               |
| (From Merge_spe_prot.py)     | **Proteic_prot_spe.treefile** |
|                              |                               |
| Proteic_prot_spe.fasta       |                               |
|                              |                               |
| (From Merge_spe_prot.py)     |                               |
+------------------------------+-------------------------------+

: Input-Output Protdist.py

## Instructions to use

### Informations_file.txt

This file is the Execution input file. You need to fill up all the information in this file in order to use the pipeline.

**!If an information is missing it may not work!**

### Separator

DO NOT USE "," AS A SEPARATOR INSIDE EACH COLONE USE ";" INSTEAD

### Example

Example for the study Orthobunyaviridae with a specific outgroup

#### A : Execution Project name

This name will be used as the folder name. It also allows you to specify which project the outgroup-specific data is associated with, where necessary.

If you wish to analyse several segments or polyproteins from the same family or species, you must enter two different project names; otherwise, all the files produced by the pipeline for both projects will be stored in the same folder.

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

#### B : Study name (if Specific OG)

This name is used to distinguish specific outgroups.

If the row does not concern a specific outgroup, enter the project name. Otherwise, enter a name for the specific outgroup (ex. name).

If you wish to use different specific outgroups for the same project,you must have one row per outgroup.

+-------------------+-------------------+
| A                 | Study name        |
+===================+===================+
| Orthobunyaviridae | Khurdun           |
+-------------------+-------------------+
| Orthobunyaviridae | Pacuvirus         |
+-------------------+-------------------+
| Orthobunyaviridae | Orthobunyaviridae |
+-------------------+-------------------+

: Study name : Where Khurdun and Pacuvirus will be used as outgroups

#### C : ID files

The name of the file that contains the virus IDs. This file must be in the same folder as the pipeline.

+-------------------+-------------------+----------------------------+
| A                 | B                 | ID files                   |
+===================+===================+============================+
| Orthobunyaviridae | Khurdun           | Khurdun_file.txt           |
+-------------------+-------------------+----------------------------+
| Orthobunyaviridae | Pacuvirus         | Pacuvirus_file.txt         |
+-------------------+-------------------+----------------------------+
| Orthobunyaviridae | Orthobunyaviridae | Orthobunyaviridae_file.txt |
+-------------------+-------------------+----------------------------+

: ID files

#### D : Outgroup Spe (Y/N)

Y : if this line concern a specific outgroup for one or multiple proteines.

N : If this line doesn't concern a specific outgroup.

+-------------------+-------------------+---------------+---------------+
| A                 | B                 | C             | Outgroup Spe  |
+===================+===================+===============+===============+
| Orthobunyaviridae | Khurdun           | IDk           | Y             |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Pacuvirus         | IDp           | Y             |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | IDn           | N             |
+-------------------+-------------------+---------------+---------------+

: Outgroup Specific : Is it a specifc outgroup?

#### E : For wich prot

Which proteins are concerned by this study/line?

If this line doesn't concern a specific outgroup, enter the names of the proteins from the studied species .

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

: For wich prot

#### F : Proteins name (in order)

Proteins name of the study name (col B)

+-------------------+-------------------+-------+--------------------------+
| A                 | B                 | E     | Proteins name (in order) |
+===================+===================+=======+==========================+
| Orthobunyaviridae | Khurdun           | nsmgn | polyprotein;gn;gc        |
+-------------------+-------------------+-------+--------------------------+
| Orthobunyaviridae | Pacuvirus         | prots | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+-------+--------------------------+
| Orthobunyaviridae | Orthobunyaviridae | prots | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+-------+--------------------------+

: Proteins name (in order)

#### G : spe prot

Identify the proteins that will be involved in the specific protein

+-------------------+-------------------+-------+---------------+
| A                 | B                 | F     | spe prot      |
+===================+===================+=======+===============+
| Orthobunyaviridae | Khurdun           | prots | nsmgn=gn      |
+-------------------+-------------------+-------+---------------+
| Orthobunyaviridae | Pacuvirus         | prots |               |
+-------------------+-------------------+-------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | prots | nsmgn=nsm+gn  |
+-------------------+-------------------+-------+---------------+

: spe prot

The protein Gn of Khurdun will be used as the outgroup.

#### H : TargetID

Enter an reference ID that will be used for BLAST alignment in case some positions are missing.

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

#### I : reference proteins position

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

#### J : Start of the proteins in the nucleic sequence

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

**K : DNAdist Threshold**

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

#### L : Organism to keep ID

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

#### M : Organism with.out a particular protein ID

+-------------------+-------------------+-----+-------------------------------------------+
| A                 | B                 | L   | Organism with.out a particular protein ID |
+===================+===================+=====+===========================================+
| Orthobunyaviridae | Khurdun           |     | OUT\*NSM=NC_078560.1                      |
+-------------------+-------------------+-----+-------------------------------------------+
| Orthobunyaviridae | Pacuvirus         |     | NONE\*NN=NN                               |
+-------------------+-------------------+-----+-------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDs | NONE\*NN=NN                               |
+-------------------+-------------------+-----+-------------------------------------------+

: Organism with.out a particular protein ID

OUT\*\<prot\>=list of IDs

IN\*\<prot\>=liste of IDs

**N : Trimal**

Do you want to use trimal for this study.

+-------------------+-------------------+-------------+--------+
| A                 | B                 | M           | Trimal |
+===================+===================+=============+========+
| Orthobunyaviridae | Khurdun           | OUT\*NSM=ID | Yes    |
+-------------------+-------------------+-------------+--------+
| Orthobunyaviridae | Pacuvirus         | NONE        | Yes    |
+-------------------+-------------------+-------------+--------+
| Orthobunyaviridae | Orthobunyaviridae | NONE        | Yes    |
+-------------------+-------------------+-------------+--------+

: Do you want to use Trimal

O : Blast

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

Please pay attention to the ***Final_pos_prot.txt*** file. The blast alignement is only an approximation

This is not a completely reliable solution. We strongly advise against using Vertex if you have very few organisations with records in the databases.

It is also strongly recommended that you check the coherence between the actual positions of the proteins and the predicted positions. The **Pos_prot.txt** file contains all the positions extracted from the databases.

If you are not satisfied with the positions, you can change them manually. After that, simply restart the pipeline from the **mark_spe.py** file.

### Region name in the database

The same proteins may have different region names from one organism to another.

If you have very few results for a protein, and you are certain that this is not normal, it may be because the protein’s name has not been recognised.

You can use the file ***prot_all\_\<Study\>.gp*** to help you identified the different region name.

Whether it differs from the name you entered in the various protein names (red rectangle).

![Gc and Gn different name](images/Diff_G1_G2.pdf)

Or that it differs from the various names of the same protein (blue rectangles)

![Polyprotein different name](images/diff_nom_polyprot.pdf)

You can add protein name equivalents to the file **protein_equivalent.txt**

such as : ***glycoproteinprecursor\<\\t\>polyprotein***

### Host/vectors caption

The names of the hosts/vectors are not all uniform and some are not recognised by the NCBI databases. To get past this, you can add equivalents for the host/vector name and the scientific name in the **species_equivalent.txt**

such as : ***mosquito\<\\t\> Culicidae***

The Images of the various orders of hosts/vectors are available on another public GitHub repository.

Feel free to add new illustrations, this GitHub repository does not yet contain all the orders.

The images were generated by AI.
