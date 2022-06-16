# rRNA-depleation-project

Welcome to [naam]

Created by Charlotte Bommelijn and Nikki rovers

This Snakemake pipeline makes a set of guide rRNA’s for Cas-9 digestion of the rRNA of a given species. It uses barrnap and bowtie in combinations with two custom python scripts to create a set of guide rna's to digest all the rRNA genes for a given species.

input a fasta file of the genomic sequence and get a list of sequences to order compatible with EnGen® sgRNA Synthesis Kit, S. pyogenes


## Installation:

#### Import the github repository containing the following files and directory’s:

- README.txt
- snakefile  
   file containing the Snakemake pipeline
- scripts/  
   directory containing all the custom python scripts for the pipeline
- genomes/  
   empty directory in with to place the genome fasta files of the species of interest, formatted as: species_name.fna
  

#### Create a conda environment with the following tools:

Barrnap version 0.9

Bowtie version 1.2.2


## Input


The whole genome sequence of the organism of interest has to placed in the genomes folder in the installation. This needs to be unzipped and the file renamed to the scientific name of the species with an underscore ( _ ) as the space and with the file extension .fna

_good examples:_

aspergillus_niger.fna  
homo_sapiens.fna  
escherichia_coli.fna


_**wrong** examples:_ 

aspergillus niger.fna  
c.auris.fasta  
mouse.fna


## Running


With the conda environment containing barrnap and bowtie active, the Snakemake pipeline can be run  
For this a name for the set needs to be given in the place of “example_set”, this name can not include a spaces.

It is good practice to test the installation with a snakemake dry run using the following command:
```
$ snakemake --cores 1 --config set=example_set --dry-run
```
Then the pipeline can be run using the following command:
```
$ snakemake  --cores 1 --config set=example_set
```
The --cores can be changed to set the number of cores you want to use. With only one the longest step is the barrnap rRNA extraction.


## Output

The pipeline outputs 2 files:
- example_set.tsv
- example_set.order.txt

The .tsv file contains in the first row a list of the species that the set is designed for. The second row shows the number of guides in the set. The following rows list the guideID, the sequence of the guide, the number of genes that it targets and a list of which genes these are. The guideID indicated the species that is was designed on but this does not mean that it only targets that specific species. At the end the file list the number of rRNA genes that no guides map to and list those sequences.

The .order.txt file contains a list of the oligo’s that need to be ordered. These do not have the pam site in the sequence and have the T7 promotor at the 5’ end and a 14nt overlap at the 3’ end making this set compatible with the [EnGen® sgRNA Synthesis Kit, S. pyogenes](https://international.neb.com/products/e3322-engen-sgrna-synthesis-kit-s-pyogenes#Product%20Information)

snakemake also creates the folowing directorys with all the files of the annalisys, these do not need to be removed or emptied between runs.

- logs/  
   output directory for the log reports and error messages for all the Snakemake rules
- guide_finder/  
   output directory with the fasta files containing all the guides on each sequence
- rRNA_fasta/  
output directory for barrnap containing the fasta sequences of the rRNA
- rRNA_gff/  
  output directory for barrnap containing the fasta sequences of the rRNA
  


## Troubleschooting

For some species barrnap is not able to identify the rRNA, the output files then are empty. The genome of this species should be removed from the genomes folder and the snakemake pipeline can be run again.

