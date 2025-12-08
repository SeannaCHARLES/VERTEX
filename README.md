---
editor_options: 
  markdown: 
    wrap: 72
---

# Script for phylogenetic studies

## Project description

Le but final de ce script est de générer des arbres phylogénétique afin
que vous puissiez les analyser et les interprétés. Ce script permet de :

### Main.sh

-   Download all the sequences in the Genbank format.

This step download all the sequences in a Genbank format concatenated
all the Genbank file into one and delete the single sequences Genbank
file.

| Input                           | Output           |
|---------------------------------|------------------|
| All_information.txt (From user) | All_sequences.gb |

: Input-Output Main.sh

### Genbank_extract.py

-   Seperate all the nucleic sequences and the proteic sequences into
    single files

    From the Genbank file the nucleic and proteic sequences are being
    extract and save in 2 folders ***nucl*** and ***prot***.

-   Create a **Metadata_files.txt** that contains different
    informations.

This step create a file that contains all the informations available in
the genbank file format and the Gp file format (the gp file format are
being download and delete after). This file contains the following
information : Study–Virus name–ID itol tree–Taxonomy–Serogroup–Protein
accesion number–Size Proteic sequence–Startproteins–Endproteines–Nucleic
accession number–Size Nucleic sequence–Strand–Strain–Geographical
region–Host or vector References

-   Create the caption for the Serogroup, Geographical region, and Host.

    The caption are created based on the Itol format.

    One color will be associated with the same values for each studies.

    The Serogroup captions are created using the **Serogroup.txt** file.
    You have to provide you own **Serogroup.txt** file since this
    information cannot be found in the genbank files format

+---------------------------------+---------------------------------------+
| Input                           | Output                                |
+:===============================:+=======================================+
| All_sequences.gb (From Main.sh) | proteic_sequence.fasta                |
|                                 |                                       |
|                                 | nucleic_sequence.fasta                |
|                                 |                                       |
|                                 | All_proteic_sequences.fasta           |
|                                 |                                       |
|                                 | All_nucleic_sequences.fasta           |
|                                 |                                       |
|                                 | All_proteic_sequences_align.aln-fasta |
|                                 |                                       |
|                                 | All_nucleic_sequences_align.aln-fasta |
|                                 |                                       |
|                                 | Serogroup_caption.txt                 |
|                                 |                                       |
|                                 | Geographical_caption.txt              |
|                                 |                                       |
|                                 | Host_caption.txt                      |
|                                 |                                       |
|                                 | Metadata.txt                          |
|                                 |                                       |
|                                 | Pos_prot.txt                          |
+---------------------------------+---------------------------------------+

: Input-Output Genbank.py

### Prot_position.py

-   This file allows the reference proteins to be aligned(using blastp)
    with proteins whose positions are missing.

This file ensures that as many protein positions as possible are
retrieved in cases where the information is not available in the Gp file
format.

| Input                          | Output             |
|--------------------------------|--------------------|
| Pos_prot.txt (From Genbank.py) | Final_pos_prot.txt |

: Input-Output Prot_position.py

### Main.sh

-   Create a Blast database of all the reference virus/organisms
    proteins

| Input                                   | Output                     |
|-----------------------------------------|----------------------------|
| Ref_prot.fasta(from Genbank_extract.py) | Proteic_prot_1-n_ref.fasta |

: Input-Output Main.sh

### DNA_dist.sh

-   This file lauch DNA dist from phylip and create a divergence
    sequence table.

This file used

| Input | Output |
|----|----|
| All_nucleic_sequences_align.aln-fasta (From Main.sh) | Divergence_nucleic_table.txt |

: Input-Output DNA_dist.sh

### Delete_dbl.py

-   This file delete de sequences that are too close using the
    **DNA_dist.sh** output table

+----------------------------------+----------------------------------+
| Input                            | Output                           |
+==================================+==================================+
| Divergence_nucleic_table.txt     | All_nucleic_sequences_sdbl.fasta |
| (From DNA_dist.sh)               |                                  |
|                                  | All_proteic_sequences_sdbl.fasta |
+----------------------------------+----------------------------------+

: Input-Output Delete_dbl.py

### Test_script_mark_spe.py

-   Seperate the different proteins and save them

    One file for each proteins.

+--------------------------------------------+------------------------+
| Input                                      | Output                 |
+============================================+========================+
| All_nucleic_sequences_sdbl.fasta (From     | Nucleic_prot_1-n.fasta |
| Delete_dbl.py)                             |                        |
|                                            | Proteic_prot_1-n.fasta |
| All_proteic_sequences_sdbl.fasta (From     |                        |
| Delete_dbl.py)                             | Nucleic_prot_OG.fasta  |
|                                            |                        |
|                                            | Proteic_prot_OG.fasta  |
+--------------------------------------------+------------------------+

: Input-Output Test_sript_mark_spe.py

### Merge_spe_prot.py

-   This file merge the specific OG species, with the specific study.

+--------------------------------------------+------------------------+
| Input                                      | Output                 |
+============================================+========================+
| Nucleic_prot_1-n.fasta (From               | Nucleic_prot_spe.fasta |
| Test_sript_mark_spe.py)                    |                        |
|                                            | Proteic_prot_spe.fasta |
| Proteic_prot_1-n.fasta (From               |                        |
| Test_sript_mark_spe.py)                    |                        |
|                                            |                        |
| Nucleic_prot_OG.fasta (From                |                        |
| Test_sript_mark_spe.py)                    |                        |
|                                            |                        |
| Proteic_prot_OG.fasta (From                |                        |
| Test_sript_mark_spe.py)                    |                        |
+--------------------------------------------+------------------------+

: Input-Output Merge_spe_prot.py

### Protdist.py

-   This file lauch Protdist from phyllip and create a sequences
    divergence table

+---------------------------------------+-----------------------------+
| Input                                 | Output                      |
+=======================================+=============================+
| Nucleic_prot_1-n.fasta(From           | Divergence_proteic_table.py |
| Test_sript_mark_spe.py)               |                             |
|                                       |                             |
| Proteic_prot_1-n.fasta(From           |                             |
| Test_sript_mark_spe.py)               |                             |
|                                       |                             |
| Nucleic_prot_spe.fasta(From           |                             |
| Merge_spe_prot.py)                    |                             |
|                                       |                             |
| Proteic_prot_spe.fasta(From           |                             |
| Merge_spe_prot.py)                    |                             |
+---------------------------------------+-----------------------------+

: Input-Output Protdist.py

### Main.sh

-   Trimal

-   IQtree2

Identified the best model (using iqtree2)

iqtree2 was used because it can provide a test to identify the best
model for the nucleic sequences and the proteic sequences. If you only
want to study the proteic sequences you may used prottest instead .

-   Create the phylogenetic tree using the best model identified (using
    iqtree2)

+-----------------------------------------+---------------------------+
| Input                                   | Output                    |
+=========================================+===========================+
| Proteic_prot_1-n.fasta(From             | Proteic_prot_1-n.treefile |
| Test_sript_mark_spe.py)                 |                           |
|                                         | Proteic_prot_spe.treefile |
| Proteic_prot_spe.fasta(From             |                           |
| Merge_spe_prot.py)                      |                           |
+-----------------------------------------+---------------------------+

: Input-Output Protdist.py

## Instructions to use

Informations_file.txt

This file is the main input file. You need to fill up all the
information in this file in order to use the package.

If an information is missing it will not work.

A : Main Study name

B : Study name (if Specific OG)

C : List of ID

D : Outgroup Spe (Y/N)

Is this study name (B) is a specific Outgroup

E : For quel prot

Quels sont les protéines concernées par cette étude

F : Proteins of this line

G : spe prot

Spécifique proteines if needed

H : TargetID

References Target ID

I : Target proteins position

J : Start of the proteins in the nucleic sequence

K : DNAdist Threshold

L : Organism to keep ID

M : Organism with(out) a particular protein ID

N : Trimal

Utiliser trimal pour cette étude?

O : Different prot name

## Instructions to adapt the code to your needs

This scripts is currently working for the peribunyaviridae (96 species)
the nairoviridae(43 species) the CCHFV strains(123 strains) and the
oropouche strains( strains) if you want to use this script or add other
organisms you need to modifiy the following points :

**The package won't work if you do not enter this information**

**If you want to add another species it is recommanded to rerun the
entire script in order to have the same caption for all your studied
species**

### If you are not using the serogroup

### If you want to use these script for other species

### If you only want to create the itol legend

### If you only want to create the Metadata files

## Pré-requis

python iqtree2 package python phyml \## Installation \## Output files

### Metadata_files.txt

This file contains the following informations(if they are available in
the NCBI genbank file format and the gb file format) on all the organims
you choose.

```         
    Virus name  ID itol tree    Taxonomy    Serogroup   Protein accesion number Size Proteique sequence Start MLD   End MLD Sart GP38   End GP38    Start Gn    End Gn  Start NSm   End NSm Start Gc    End Gc  Nucleique accession number  Size Nucelique sequence Strand  Strain  Geographical region Host    Host source References
```
