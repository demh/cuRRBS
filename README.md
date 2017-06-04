# cuRRBS

Created by Daniel E. Martin-Herranz, Antonio J.M. Ribeiro and Thomas M. Stubbs.

Copyright (C) 2016,2017 D.E. Martin-Herranz, A.J.M. Ribeiro, T.M. Stubbs.

## What is cuRRBS?

**cuRRBS** (customised Reduced Representation Bisulfite Sequencing) is an easy-to-use computational method that predicts what is the optimal combination of *restriction enzymes* and *size range* that should be used in order to enrich for any given set of sites of interest in any genome. In other words, by modifying two steps of the original protocol, cuRRBS generalises [RRBS](http://www.nature.com/nprot/journal/v6/n4/full/nprot.2010.190.html). This allows to create 'personalised' reduced representations (i.e. subsets) of the genome, that make DNA methylation next-generation sequencing experiments more cost-effective.

## Getting started

Before being able to run cuRRBS in your computer, you will need to download and install the following software:

1. [Python 2.7](https://www.python.org/downloads/). The programming language in which cuRRBS is written.
2. [NumPy](https://scipy.org/install.html#individual-packages). A Python module which contains many useful functions used by cuRRBS.
3. [Terminal](https://en.wikipedia.org/wiki/Comparison_of_terminal_emulators). You will need it to run cuRRBS commands. 

Afterwards, you can install cuRRBS in your computer by [direct download](https://github.com/demh/cuRRBS/archive/master.zip) or cloning it from this GitHub repository:

```
git clone https://github.com/demh/cuRRBS.git
```

## Preparing cuRRBS input files

cuRRBS requires the following files as input:

**1. Isoschizomer annotation file.** 

Restriction enzymes that cut the genome in the same way (i.e. generate the same fragment length distribution) are grouped in isoschizomer families. Furthermore, this file also contains information regarding the methylation-sensitivity (in a concrete genomic context) of the commercially-available restriction enzymes. These files can be found in utils/isoschizomers_*_annotation.csv but the user can also create his own file. By default, cuRRBS uses isoschizomers_CpG_annotation.csv, which should be used when considering only methylation-sensitivity in CpG contexts (i.e. for most vertebrate genomes). 

**2. Pre-computed files.** 

They contain information for the *in silico* digestions of a genome with certain restriction enzymes. There are two ways to obtain these files:

* [Download them](http://www.ebi.ac.uk/~dem44/cuRRBS_pre_computed_files/) (recommended option for human, mouse and *Arabidopsis thaliana* genomes). Only the enzymes which belong to an isoschizomer family with CpG methylation-insensitive members (see utils/isoschizomers_CpG_annotation.csv) are available. They need to be extracted using the following command in the terminal:

   ```
   tar xvjf pre_computed_files_for_genome_of_interest.tar.bz2 /path/to/precomputed/files/folder/
   ```

* Generated by the user (recommended only when using a novel genome). The script src/pre_compute_digestions.py can be used for these purposes:

   ```   
   python pre_compute_digestions.py enzymes_to_pre_compute.txt genome.fa working_directory
   ``` 

**3. Enzymes to check.** 

This file contains the isoschizomer families that should be considered by cuRRBS (i.e. those that contain at least one methylation-insensitive member). When working with CpG methylation-insensitive enzymes, the file utils/enzymes_to_check_CpG.txt should be used.

**4. Sites annotation file.** 

This file stores the information regarding the sites of interest (i.e. genomic coordinates) that the user wants to enrich for with the new customised protocol. This file needs to be created from scratch by the user, since it will be different depending on the sites that need to be targeted. The file is provided with a CSV (comma-separated values) format, with each row containing the information for a site of interest and the following columns:  

* Site_ID: unique ID for the current site of interest. It can not contain the '_' or the ',' characters.

* Chr: chromosome where the current site of interest is located. The same nomenclature should be used as in the case of the pre-computed files.

* Coordinate: 1-based genomic coordinate where the current site of interest is located. Please make sure that the coordinates are based on the same genome assembly as the one used to generate the pre-computed files.

* Weight: in case there is a preference in recovering certain sites of interest, their weights should be higher. The weights are always positive and are used afterwards to calculate the *Score* (see 'Interpreting cuRRBS output').

The first line of the sites annotation file must be a header containing the column names. Some examples of these type of files for different biological systems can be found in the examples/ folder.


   
