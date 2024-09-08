# slice_multi_gbk
Slices/subsections genbank files by either gene names, locus tags or base pair range. This handles multi-genbank files which have multiple contigs. Outputs sub-setted genbank files as well as multifasta CDS and amino acid files.

## Requirements
- biopython >1.84

## Installation
Easiest to create a mamba/conda environment
```bash
mamba create -y -n slice_gbk_env
mamba install -y -n slice_gbk_env conda-forge::biopython=1.84
```
Test script
```bash
mamba activate slice_gbk_env
python slice_multi_gbk.py -h
```

## Quick start
```bash
# Extract the region between two genes including these genes
python slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output # using gene names
python slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output -p # using gene names, also outputting multifasta protein and CDS files with sliced genbank
python slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output # using locus tags

# Extract 2 genes up- and down-stream of your specified two genes
python slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n 2 -o gene_sliced_output # using gene names
python slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n 2 -o locus_tag_sliced_output # using locus tags

# Extract a single gene
python slice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output # using gene names
python slice_multi_gbk.py -i infile.gbk -lt locus_tag1: -o locus_tag_sliced_output # using locus tags

# Use a wildcard to extract gene/s
python slice_multi_gbk.py -i infile.gbk -g gene*: -o wildcard_gene_extract_output # using gene names
python slice_multi_gbk.py -i infile.gbk -lt locus_tag*: -o wildcard_gene_extract_output # using locus tags

# Extract a single gene and 3 genes up- and down-stream of the gene
python slice_multi_gbk.py -i infile.gbk -g gene1: -n 3 -o gene_extract_output # using gene names
python slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n 3 -o gene_extract_output # using locus tags

# Extract a base pair range from contig. Contig name required to handle multi-genbank files
python slice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output
```

## Usage and options

| Command   | Decsription  |   Example    |
|  --- | ---   |    -------   |
| `-i`/`--infile`  | Input genbank file. Required  |     |
| `-g`/`--genes`  | The one/two gene names to slice between (inclusive). Handles wildcards (*). | Single gene: `-g hemF:`. Between two genes: `-g hemF:maeB` |
| `-lt`/`--locus_tags`  | One/two locus tags of gene to slice between (inclusive). Handles wildcards (*) | Single gene: `-lt DKHHJM_01915:`. Between two genes: `-lt DKHHJM_01915:DKHHJM_01920`    |
| `-n`/`--num_genes`  | Number of genes to slice before and after the selected genes. Will take maximum number of genes up/downstream if at end of a contig. Used with `-g` or `-lt`. Default=0  | Two genes up and downstream: `-n 2`  |
| `-r`/`--range`  | Two base pair coordinates (range) to slice between. Contig name or locus must be provided. Used instead of `-g` or `-lt`. Optional | Slicing contig_1 between bp 612899-630276: `-r 612899:630276:contig_1`    |
| `-o`/`--outfile`  | Output directory and filename prefix. Required  | Name output E_coli: `-o E_coli`  |
| `-p`/`--protein`  | Produce multifasta files for nucleotide (CDS) and amino acid sequences contained within sliced genbanks. Optional  | Get multifasta: `-p`   |
| `-c`/`--case_insensitive`  | Turn **OFF** case sensitivity for gene and locus tag matching. Optional  | Turn off case sensitivity: `-c`    |
| `-v`/`--version`  | Get version of script. Run alone  | Get version: `-v`   |

## Authors
- Ben Vezina
   - Scholar: https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en&oi=ao
   - ORCID: https://orcid.org/0000-0003-4224-2537
- Abhinaba Ray
- Adapted from https://gist.github.com/jrjhealey/2df3c65c7a70cbca4862e94620e4d7b2
