---
title: "VERTEX README"
author: "CHARLES Seanna"
date: "2026-09-18"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# VERTEX

![Workflow of the vertex](images/Workflow_v6.pdf)

## Project description

The purpose of this pipeline is to automaticaly separate protein or proteic domain of sequences, and finally generate phylogenetic trees with the separate region.

This pipeline use the NCBI database and is therefore dependent not only on the information provided but also on the uniformity of that information.

It will create caption that you can directly use on the website Itol. It will also allow you to create other files that may be useful for your research.

### Execution.sh

- Download all the sequences in the Genbank format.

This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.

+---------------------+------------------+
| Input               | Output           |
+:===================:+:================:+
| all_information.txt | all_sequences.gb |
|                     |                  |
| (From user)         |                  |
+---------------------+------------------+

: Input-Output Execution.sh

### Genbank_extract.py

- Separate all the nucleic sequences and the proteic sequences into single files

From the Genbank file the nucleic and proteic sequences are being extract and save in 2 folders nucl and prot.

- Create a **Metadata_files.txt** that contains different information available in the genbank file format and the Gp file format.

This file contains the following information : Study–Virus name–ID itol tree–Taxonomy–Informations–Protein accession number–Size Proteic sequence–Startproteins–Endproteines–Nucleic accession number–Size Nucleic sequence–Strand–Strain–Geographical region–Host or vector References

The gp file format are being download and concatenate in the ***prot_all\_\<Study_name\>.gp*** file).

- Create the caption for the Informations, Geographical region, and Host.

The caption are created based on the Itol format.

One color will be associated with the same values for all the studies.

The Informations captions are created using the **personal_informations.txt** file. You can put in that file all the information you want to use to caracterise you study (for example serogroup or precise taxonomy)

+---------------------+---------------------------------------+
| Input               | Output                                |
+:===================:+:=====================================:+
| All_sequences.gb    | proteic_sequence.fasta                |
|                     |                                       |
| (From Execution.sh) | nucleic_sequence.fasta                |
|                     |                                       |
|                     | All_proteic_sequences.fasta           |
|                     |                                       |
|                     | All_nucleic_sequences.fasta           |
|                     |                                       |
|                     | All_proteic_sequences_align.aln-fasta |
|                     |                                       |
|                     | All_nucleic_sequences_align.aln-fasta |
|                     |                                       |
|                     | All_proteic_sequences.gp              |
|                     |                                       |
|                     | **Informations_caption.txt**          |
|                     |                                       |
|                     | **Geographical_caption.txt**          |
|                     |                                       |
|                     | **Host_caption.txt**                  |
|                     |                                       |
|                     | **Metadata.txt**                      |
|                     |                                       |
|                     | Pos_prot.txt                          |
+---------------------+---------------------------------------+

: Input-Output Genbank.py

### Prot_position.py (optional)

- This file use blastp to complete the missing position within a study.

This file ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file format.

+-------------------+------------------------+
| Input             | Output                 |
+:=================:+:======================:+
| Pos_prot.txt      | **Final_pos_prot.txt** |
|                   |                        |
| (From Genbank.py) |                        |
+-------------------+------------------------+

: Input-Output Prot_position.py

### Execution.sh

- Create a Blast database of all the reference virus/organisms proteins.

+---------------------------+----------------------------+
| Input                     | Output                     |
+:=========================:+:==========================:+
| Ref_prot.fasta            | Proteic_prot_1-n_ref.fasta |
|                           |                            |
| (From Genbank_extract.py) |                            |
+---------------------------+----------------------------+

: Input-Output Execution.sh

### DNA_dist.sh

- This file lauch DNA dist from phylip and create a divergence sequence table.

This file used

+-----------------------------------------+----------------------------------+
| Input                                   | Output                           |
+:=======================================:+:================================:+
| All \_nucleic_sequences_align.aln-fasta | **Divergence_nucleic_table.txt** |
|                                         |                                  |
| (From Execution.sh)                     |                                  |
+-----------------------------------------+----------------------------------+

: Input-Output DNA_dist.sh

### Delete_dbl.py

- This file delete de sequences that are judged too close using the **DNA_dist.sh** output table and the define threshold define in the **All_information.txt**.

+------------------------------+----------------------------------+
| Input                        | Output                           |
+:============================:+:================================:+
| Divergence_nucleic_table.txt | All_nucleic_sequences_sdbl.fasta |
|                              |                                  |
| (From DNA_dist.sh)           | All_proteic_sequences_sdbl.fasta |
+------------------------------+----------------------------------+

: Input-Output Delete_dbl.py

### Test_script_mark_spe.py

- Separate the different proteins and save them

One file for each proteins.

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

- This file lauch Protdist from phyllip and create a sequences divergence table for each protein.

+-------------------------------+---------------------------------+
| Input                         | Output                          |
+:=============================:+:===============================:+
| Nucleic_prot_1-n.fasta        | **Divergence_proteic_table.py** |
|                               |                                 |
| (From Test_sript_mark_spe.py) |                                 |
|                               |                                 |
| Proteic_prot_1-n.fasta        |                                 |
|                               |                                 |
| (From Test_sript_mark_spe.py) |                                 |
|                               |                                 |
| Nucleic_prot_spe.fasta        |                                 |
|                               |                                 |
| (From Merge_spe_prot.py)      |                                 |
|                               |                                 |
| Proteic_prot_spe.fasta        |                                 |
|                               |                                 |
| (From Merge_spe_prot.py)      |                                 |
+-------------------------------+---------------------------------+

: Input-Output Protdist.py

### Execution.sh

- Trimal

- IQtree2

Identified the best model (using iqtree2)

iqtree2 was used because it can provide a test to identify the best model for the nucleic sequences(not present in this pipeline) and the proteic sequences. If you only want to study the proteic sequences you may used prottest instead.

- Create the phylogenetic tree using the best model identified (using iqtree2)

+-------------------------------+-------------------------------+
| Input                         | Output                        |
+:=============================:+:=============================:+
| Proteic_prot_1-n.fasta        | **Proteic_prot_1-n.treefile** |
|                               |                               |
| (From Test_sript_mark_spe.py) | **Proteic_prot_spe.treefile** |
|                               |                               |
| Proteic_prot_spe.fasta        |                               |
|                               |                               |
| (From Merge_spe_prot.py)      |                               |
+-------------------------------+-------------------------------+

: Input-Output Protdist.py

## Instructions to use

Informations_file.txt

This file is the Execution input file. You need to fill up all the information in this file in order to use the pipeline.

If an information is missing it will not work!

**DO NOT USE "," AS A SEPARATOR INSIDE EACH COLONE USE ";" INSTEAD**

Example for the study Orthobunyaviridae with a specific outgroup

A : Execution Project name

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

B : Study name (if Specific OG)

This name is used to distinguish specific outgroups.

If the row does not concern a specific outgroup, enter the project name. Otherwise, enter a name for the specific outgroup (its name, for example).

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

C : ID files

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

D : Outgroup Spe (Y/N)

Y : if your project have a specific outgroup for one or multiple proteines.

N : If your project doesn't have a specific outgroup.

+-------------------+-------------------+-------------+--------------+
| A                 | B                 | C           | Outgroup Spe |
+===================+===================+=============+==============+
| Orthobunyaviridae | Khurdun           | IDk         | Y            |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Pacuvirus         | IDp         | Y            |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Orthobunyaviridae | IDn         | N            |
+-------------------+-------------------+-------------+--------------+

: Outgroup Specific : Is it a specifc outgroup?

E : For wich prot

Which proteins are covered by this study/line?

If this line refers to the project line, enter the names of the proteins from the species studied.

If this line refers to specific outgroups, specify for which proteins this outgroup will be used

+-------------------+-------------------+-------------+-----------------------+
| A                 | B                 | D           | For wich prot         |
+===================+===================+=============+=======================+
| Orthobunyaviridae | Khurdun           | Y           | nsmgn                 |
+-------------------+-------------------+-------------+-----------------------+
| Orthobunyaviridae | Pacuvirus         | Y           | polyprotein;gn;nsm;gc |
+-------------------+-------------------+-------------+-----------------------+
| Orthobunyaviridae | Orthobunyaviridae | N           | polyprotein;gn;nsm;gc |
+-------------------+-------------------+-------------+-----------------------+

: For wich prot

F : Proteins name (in order)

Proteins name of the study name (col B)

+-------------------+-------------------+---------------+--------------------------+
| A                 | B                 | E             | Proteins name (in order) |
+===================+===================+===============+==========================+
| Orthobunyaviridae | Khurdun           | nsmgn         | polyprotein;gn;gc        |
+-------------------+-------------------+---------------+--------------------------+
| Orthobunyaviridae | Pacuvirus         | prots         | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+---------------+--------------------------+
| Orthobunyaviridae | Orthobunyaviridae | prots         | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+---------------+--------------------------+

: Proteins name (in order)

G : spe prot

Identify the proteins that will be involved in the formation of the specific protein

+-------------------+-------------------+-------------+--------------+
| A                 | B                 | F           | spe prot     |
+===================+===================+=============+==============+
| Orthobunyaviridae | Khurdun           | prots       | nsmgn=gn     |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Pacuvirus         | prots       |              |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Orthobunyaviridae | prots       | nsmgn=nsm+gn |
+-------------------+-------------------+-------------+--------------+

: spe prot

H : TargetID

Enter an organisation ID to be used as a reference for the BLAST alignment.

+-------------------+-------------------+---------------+---------------+
| A                 | B                 | G             | RefID         |
+===================+===================+===============+===============+
| Orthobunyaviridae | Khurdun           | nsmgn=gn      |               |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Pacuvirus         |               |               |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | nsmgn=nsm+gn  | NC_001926.1   |
+-------------------+-------------------+---------------+---------------+

: RefID

I : reference proteins position

The positions of the protein of the reference or the specific OG.

+-------------------+-------------------+-------------+---------------------------------+
| A                 | B                 | H           | reference /OG proteins position |
+===================+===================+=============+=================================+
| Orthobunyaviridae | Khurdun           |             | 23;292;526;944                  |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Pacuvirus         |             | 1;1458;23;305;308;476;598;1386  |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDref       | 22;302;303;476;509;1372         |
+-------------------+-------------------+-------------+---------------------------------+

: RefID

J : Start of the proteins in the nucleic sequence

+-------------------+-------------------+------------+-----------------------------------------------+
| A                 | B                 | I          | Start of the proteins in the nucleic sequence |
+===================+===================+============+===============================================+
| Orthobunyaviridae | Khurdun           | pos        | 34                                            |
+-------------------+-------------------+------------+-----------------------------------------------+
| Orthobunyaviridae | Pacuvirus         | pos        | 53                                            |
+-------------------+-------------------+------------+-----------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | pos        | 57                                            |
+-------------------+-------------------+------------+-----------------------------------------------+

: Start of the proteins in the nucleic sequence

K : DNAdist Threshold

+-------------------+-------------------+-------------+-------------------+
| A                 | B                 | J           | DNAdist Threshold |
+===================+===================+=============+===================+
| Orthobunyaviridae | Khurdun           | 34          | 0.01              |
+-------------------+-------------------+-------------+-------------------+
| Orthobunyaviridae | Pacuvirus         | 53          | 0.05              |
+-------------------+-------------------+-------------+-------------------+
| Orthobunyaviridae | Orthobunyaviridae | 57          | 0.01              |
+-------------------+-------------------+-------------+-------------------+

: DNAdist Threshold

L : Organism to keep ID

organism you don't want to be deleted.

+-------------------+-------------------+-------------+---------------------------------+
| A                 | B                 | k           | Organism to keep ID             |
+===================+===================+=============+=================================+
| Orthobunyaviridae | Khurdun           | 0.01        |                                 |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Pacuvirus         | 0.05        |                                 |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | 0.01        | MH484340.1;NC_005775;EF990927.1 |
+-------------------+-------------------+-------------+---------------------------------+

: ID Organism to keep

M : Organism with.out a particular protein ID

+-------------------+-------------------+------------+-------------------------------------------+
| A                 | B                 | L          | Organism with.out a particular protein ID |
+===================+===================+============+===========================================+
| Orthobunyaviridae | Khurdun           |            | OUT\*NSM=NC_078560.1                      |
+-------------------+-------------------+------------+-------------------------------------------+
| Orthobunyaviridae | Pacuvirus         |            | NONE\*NN=NN                               |
+-------------------+-------------------+------------+-------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDs        | NONE\*NN=NN                               |
+-------------------+-------------------+------------+-------------------------------------------+

: Organism with.out a particular protein ID

N : Trimal

Do you want to use trimal for this study.

+-------------------+-------------------+-------------+-------------+
| A                 | B                 | M           | Trimal      |
+===================+===================+=============+=============+
| Orthobunyaviridae | Khurdun           | OUT\*NSM=ID | Yes         |
+-------------------+-------------------+-------------+-------------+
| Orthobunyaviridae | Pacuvirus         | NONE        | Yes         |
+-------------------+-------------------+-------------+-------------+
| Orthobunyaviridae | Orthobunyaviridae | NONE        | Yes         |
+-------------------+-------------------+-------------+-------------+

: Organism with.out a particular protein ID

O : Blast?

Do you want to use blast in case of missing positions

+-------------------+-------------------+----------+----------+
| A                 | B                 | N        | Blast    |
+===================+===================+==========+==========+
| Orthobunyaviridae | Khurdun           | Y        | Y        |
+-------------------+-------------------+----------+----------+
| Orthobunyaviridae | Pacuvirus         | Y        | Y        |
+-------------------+-------------------+----------+----------+
| Orthobunyaviridae | Orthobunyaviridae | Y        | Y        |
+-------------------+-------------------+----------+----------+

: Organism with.out a particular protein ID

### Metadata_files.txt

This file contains the following information(if they are available in the NCBI genbank file format and the gb file format) on all the organisms you choose.

``` txt
    Virus name  ID itol tree    Taxonomy    Informations   Protein accesion number Size
    Proteique sequence Start MLD   End MLD Sart GP38   End GP38    Start Gn    End Gn 
    Start NSm   End NSm Start Gc    End Gc  Nucleique accession number  Size Nucelique
    sequence Strand  Strain  Geographical region Host    Host source References
```

## Additional Instructions and warnings

### Proteins positions

Please pay attention to the ***Final_pos_prot.txt*** file. The blast alignement is only an approximation ce n'est pas une solution complétement fiable. Il est fortement déconseillé d'utiliser vertex si vous n'avez que très peu d'organisme avec des position dans les bases de données.

Il est également fortement recommandé de vérifier la cohérence entre les positions réelle de protéine et les positions prédite. Le fichier ***Pos_prot.txt***contient toutes les positions extraites des base de données.

Si les positions ne vous conviennent pas vous pouvez les modifier manuellement. Après quoi il vous suffit de relancer le pipeline à partir du ***script_mark_spe.py.***

### Region name in the database

A same proteins may have different region name d'un organismes à l'autre.

Si vous avez très peu de positions pour une protéines cela est peut-être dû à une "non reconnaissance du nom de la protéine.

You can use the file ***prot_all\_\<Study\>.gp*** to help you identified the different region name.

Que ce soit différent par rapport au nom que vous avez saisie dans les différents nom des proteines red rectangles.

![Gc and Gn different name](images/Diff_G1_G2.pdf)

ou que soit différent par rapport aux différents noms de la même protéine. blue rectangles

![Polyprotein different name](images/diff_nom_polyprot.pdf)

### Host/vectors caption

Images of the various orders of hosts are available on another public GitHub repository. Feel free to add new illustrations, this GitHub repository does not yet contain all the orders. The images were generated by AI.

## Pré-requis

python iqtree2 pipeline python phyml \## Installation \## Output files
