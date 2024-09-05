#!/usr/bin/env python3
# Slices/subsections genbank files by either gene names or base pair range. This handles multi-genbank files which have multiple contigs

# Usage:
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output.gbk # This will extract the region between two genes including these genes
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output.gbk # This will extract the region between two locus tags if gene names not available
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n 2 -o gene_sliced_output.gbk # This will extract 2 genes up- and down-stream of your specified two genes
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n 2 -o locus_tag_sliced_output.gbk
#   slice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output.gbk # This will extract a single gene
#   slice_multi_gbk.py -i infile.gbk -g gene1: -n 3 -o gene_extract_output.gbk # This will extract a single gene and 3 genes up- and down-stream of the gene
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n 3 -o gene_extract_output.gbk # This will extract a single locus tag and 3 genes up- and down-stream of this
#   slice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output.gbk # This will slice via base pair range and locus/contig. This is required to handle multi-genbank files

# Inputs:
#   1) Genbank file [.gbk, .gbff]
#   2) One of:
#       - Two gene names
#       - Two base pair range and locus

# Requirements:
#   - biopython (any)

# Author: Ben Vezina, Abhinaba Ray
#   - Scholar: https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en&oi=ao
#   - ORCID: https://orcid.org/0000-0003-4224-2537

# Citation: https://gist.github.com/bananabenana/20ff257f237d5a6e6f449fd7066577a1
# Adapted from https://gist.github.com/jrjhealey/2df3c65c7a70cbca4862e94620e4d7b2

from Bio import SeqIO
import sys, argparse
from argparse import RawTextHelpFormatter


def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description='''
description:
  Subset genbanks between 1 gene, 2 genes or base pair ranges.\n
usage:

  To slice/extract by gene names:
  \tslice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output.gbk
  
  To slice by locus tag if gene name is not available:
  \tslice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output.gbk

  To slice/extract by gene name +- a number of genes up- and down-stream:
  \tslice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n [integer] -o gene_sliced_output.gbk

  To slice/extract by locus tag +- a number of genes up- and down-stream:
  \tslice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n [integer] -o locus_tag_sliced_output.gbk

  To extract a gene by its name:
  \tslice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output.gbk

  To extract a gene by its locus tag:
  \tslice_multi_gbk.py -i infile.gbk -lt locus_tag1: -o locus_tag_extract_output.gbk

  To extract a gene by its name +- a number of genes up- and down-stream:
  \tslice_multi_gbk.py -i infile.gbk -g gene1: -n [integer] -o gene_extract_output.gbk

  To extract a gene by its locus tag +- a number of genes up- and down-stream:
  \tslice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n [integer] -o locus_tag1_extract_output.gbk

  To slice by base pair range coordinates:
  \tslice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output.gbk''', formatter_class=RawTextHelpFormatter)

        parser.add_argument('-i', '--infile', action='store',
                            help='Input genbank [.gbk, .gbff] file to slice or subsection. Can be a multi-genbank file with multiple contigs')
        parser.add_argument('-g', '--genes', action='store',
                            help='''The two gene names to slice between. Example: -g hemF:maeB
A single gene can be extracted using one gene name. Example: -g hemF:''')
        parser.add_argument('-n', '--num_genes', help='''Number of genes to slice before and after the selected genes.
        Can be used with -g or -lt only. Will take maximum number of genes up/downstream if at end of a contig. Optional. Example: -n 2''', type=int, default=0)
        parser.add_argument('-lt', '--locus_tags', help='''The locus tag of the gene to slice if gene name is not available. Example: -lt DKHHJM_01915: ''')
        parser.add_argument('-r', '--range', action='store',
                            help='The two base pair coordinates (range) to slice between. Contig name or locus must be provided. Example: 612899:630276:contig_1')
        parser.add_argument("-o", "--outfile", required=True, help="Output GenBank filename")
        parser.add_argument('-v', '--version', action='version', version='%(prog)s version 2.0.0')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")

    return parser.parse_args()


def slice_gene(record, genes, num_genes):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    try:
        # Find indices of the target genes
        target_indices = [idx for idx, feat in enumerate(loci) if 'gene' in feat.qualifiers and feat.qualifiers['gene'][0] in genes]
        
        if not target_indices:
            # print(f"No matching genes found for {genes}")
            return None

        min_target_idx = min(target_indices)
        max_target_idx = max(target_indices)

        # Calculate start index, ensuring it doesn't go below 0
        start_idx = max(0, min_target_idx - num_genes)
        
        # Calculate end index, ensuring it doesn't exceed the list length
        end_idx = min(len(loci) - 1, max_target_idx + num_genes)

        start_locus = loci[start_idx]
        end_locus = loci[end_idx]

        start = int(start_locus.location.start)
        end = int(end_locus.location.end)

        subrecord = record[start:end]
        
        actual_upstream = min_target_idx - start_idx
        actual_downstream = end_idx - max_target_idx
        
        print(f"Extracted {genes[0]}:{genes[-1]} with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id}")
        return subrecord
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None    

def slice_locustag(record, locus_tags, num_genes):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    try:
        # Find indices of the target locus tags
        target_indices = [idx for idx, feat in enumerate(loci) if 'locus_tag' in feat.qualifiers and feat.qualifiers['locus_tag'][0] in locus_tags]
        
        if not target_indices:
            # print(f"No matching locus tags found for {locus_tags}")
            return None

        min_target_idx = min(target_indices)
        max_target_idx = max(target_indices)

        # Calculate start index, ensuring it doesn't go below 0
        start_idx = max(0, min_target_idx - num_genes)
        
        # Calculate end index, ensuring it doesn't exceed the list length
        end_idx = min(len(loci) - 1, max_target_idx + num_genes)

        start_locus = loci[start_idx]
        end_locus = loci[end_idx]

        start = int(start_locus.location.start)
        end = int(end_locus.location.end)

        subrecord = record[start:end]
        
        actual_upstream = min_target_idx - start_idx
        actual_downstream = end_idx - max_target_idx
        
        print(f"Extracted {locus_tags[0]}:{locus_tags[-1]} with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id}")
        return subrecord
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None
    
def slice_range(record, start, end):
    subrecord = record[int(start):int(end)]
    print(f"Extracted {start}:{end} from {record.id}")
    return subrecord

def main():
    args = get_args()

    # Check that options are not used simultaneously
    if sum([bool(args.genes), bool(args.range), bool(args.locus_tags)]) > 1:
        sys.stderr.write("Error: -lt, -g and -r options cannot be used simultaneously.\n")
        exit(1)

    modified_records = []
    records_processed = 0
    records_with_match = 0

    records = SeqIO.parse(args.infile, 'genbank')

    for record in records:
        records_processed += 1
        sys.stderr.write(f'Searching {record.id}.\n')
        
        if args.genes:
            subrecord = slice_gene(record, args.genes.split(':'), args.num_genes)
            if subrecord:
                modified_records.append(subrecord)
                records_with_match += 1
        
        elif args.locus_tags:
            subrecord = slice_locustag(record, args.locus_tags.split(':'), args.num_genes)
            if subrecord:
                modified_records.append(subrecord)
                records_with_match += 1
        
        elif args.range:
            range_parts = args.range.split(':')
            if len(range_parts) != 3:
                sys.stderr.write("Error: When using -r option, provide start:end:contig.\n")
                exit(1)
            start, end, locus = range_parts[0], range_parts[1], range_parts[2]
            if record.name == locus:
                subrecord = slice_range(record, start, end)
                if subrecord:
                    modified_records.append(subrecord)
                    records_with_match += 1

    # Output summary
    if records_with_match > 0:
        print(f"Found matches in {records_with_match} out of {records_processed} records.")
        # Write modified records to output file
        with open(args.outfile, "w") as output_handle:
            SeqIO.write(modified_records, output_handle, "genbank")
        print(f"Written {len(modified_records)} records to {args.outfile}")
    else:
        if args.genes:
            print(f"No matching genes found for {args.genes} in any of the {records_processed} records processed.")
        elif args.locus_tags:
            print(f"No matching locus tags found for {args.locus_tags} in any of the {records_processed} records processed.")
        elif args.range:
            print(f"No matching range found for {args.range} in any of the {records_processed} records processed.")
        else:
            print(f"No matches found in any of the {records_processed} records processed.")

if __name__ == '__main__':
    main()