# slice_multi_gbk
Slices/subsections genbank files by either gene names, locus tags or base pair range. This handles multi-genbank files which have multiple contigs.

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
slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output.gbk # using gene names
slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output.gbk # using locus tags

# Extract 2 genes up- and down-stream of your specified two genes
slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n 2 -o gene_sliced_output.gbk # using gene names
slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n 2 -o locus_tag_sliced_output.gbk # using locus tags

# Extract a single gene
slice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output.gbk # using gene names
slice_multi_gbk.py -i infile.gbk -lt locus_tag1: -o locus_tag_sliced_output.gbk # using locus tags

# Use a wildcard to extract gene/s
slice_multi_gbk.py -i infile.gbk -g gene*: -o wildcard_gene_extract_output

# Extract a single gene and 3 genes up- and down-stream of the gene
slice_multi_gbk.py -i infile.gbk -g gene1: -n 3 -o gene_extract_output.gbk # using gene names
slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n 3 -o gene_extract_output.gbk # using locus tags

# Extract a base pair range from contig. Contig name required to handle multi-genbank files
slice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output.gbk
```

## Authors
- Ben Vezina
   - Scholar: https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en&oi=ao
   - ORCID: https://orcid.org/0000-0003-4224-2537
- Abhinaba Ray
- Adapted from https://gist.github.com/jrjhealey/2df3c65c7a70cbca4862e94620e4d7b2
