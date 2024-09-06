#!/usr/bin/env python3
# Slices/subsections genbank files by either gene names or base pair range. This handles multi-genbank files which have multiple contigs

# Usage:
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output.gbk # This will extract the region between two genes including these genes
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output.gbk # This will extract the region between two locus tags if gene names not available
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n 2 -o gene_sliced_output.gbk # This will extract 2 genes up- and down-stream of your specified two genes
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n 2 -o locus_tag_sliced_output.gbk
#   slice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output.gbk # This will extract a single gene
#   slice_multi_gbk.py -i infile.gbk -g gene*: -o gene_extract_output.gbk # This will extract any genes matching this wildcard
#   slice_multi_gbk.py -i infile.gbk -g gene1: -n 3 -o gene_extract_output.gbk # This will extract a single gene and 3 genes up- and down-stream of the gene
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n 3 -o gene_extract_output.gbk # This will extract a single locus tag and 3 genes up- and down-stream of this
#   slice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output.gbk # This will slice via base pair range and locus/contig. This is required to handle multi-genbank files

# Inputs:
#   1) Genbank file [.gbk, .gbff, .gb]
#   2) One of:
#       - Two gene names
#       - Two base pair range and locus

# Requirements:
#   - biopython (1.84)

# Author: Ben Vezina, Abhinaba Ray
#   - Scholar: https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en&oi=ao
#   - ORCID: https://orcid.org/0000-0003-4224-2537

# Citation: https://gist.github.com/bananabenana/20ff257f237d5a6e6f449fd7066577a1
# Adapted from https://gist.github.com/jrjhealey/2df3c65c7a70cbca4862e94620e4d7b2

from Bio import SeqIO
import sys, argparse
import os
from argparse import RawTextHelpFormatter
import fnmatch

def get_args():
    """Parse command line arguments"""
    try:
        parser = argparse.ArgumentParser(
            description='Subset genbanks between 1 gene, 2 genes or base pair ranges.\n'
                        'usage:\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n [integer] -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n [integer] -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1: -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene*: -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1: -n [integer] -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n [integer] -o prefix\n'
                        '  slice_multi_gbk.py -i infile.gbk -r start:end:contig -o prefix\n',
            formatter_class=RawTextHelpFormatter)

        parser.add_argument('-i', '--infile', action='store', help='Input genbank [.gbk, .gbff] file to slice or subsection. Can be a multi-genbank file with multiple contigs')
        parser.add_argument('-g', '--genes', action='store', help='The two gene names to slice between. Handles wildcards (*). Example: -g hemF:maeB. A single gene can be extracted using one gene name. Example: -g hemF:')
        parser.add_argument('-n', '--num_genes', help='Number of genes to slice before and after the selected genes. Can be used with -g or -lt only. Will take maximum number of genes up/downstream if at end of a contig. Optional. Example: -n 2', type=int, default=0)
        parser.add_argument('-lt', '--locus_tags', help='The locus tag of the gene to slice if gene name is not available. Handles wildcards (*). Example: -lt DKHHJM_01915:')
        parser.add_argument('-r', '--range', action='store', help='The two base pair coordinates (range) to slice between. Contig name or locus must be provided. Example: 612899:630276:contig_1')
        parser.add_argument('-o', '--outfile', required=True, help='Output directory and filename prefix')
        parser.add_argument('-c', '--case_insensitive', action='store_true', help='Turn off case sensitivity for gene and locus tag matching.')
        parser.add_argument('-v', '--version', action='version', version='%(prog)s version 2.1.0')

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write('An exception occurred with argument parsing. Check your provided options.')

    return parser.parse_args()

def gene_matches(gene_name, pattern, case_insensitive):
    """Check if a gene name matches the pattern, respecting case sensitivity."""
    if case_insensitive:
        return fnmatch.fnmatch(gene_name.lower(), pattern.lower())
    else:
        return fnmatch.fnmatch(gene_name, pattern)

def locus_tag_matches(locus_tag, pattern, case_insensitive):
    """Check if a locus tag matches the pattern, respecting case sensitivity."""
    if case_insensitive:
        return fnmatch.fnmatch(locus_tag.lower(), pattern.lower())
    else:
        return fnmatch.fnmatch(locus_tag, pattern)

def slice_gene(record, genes, num_genes, case_insensitive, output_prefix):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    
    try:
        target_indices = [idx for idx, feat in enumerate(loci) if 'gene' in feat.qualifiers and any(gene_matches(feat.qualifiers['gene'][0], gene, case_insensitive) for gene in genes)]
        
        if not target_indices:
            return 0

        extracted_count = 0

        for idx in range(len(target_indices)):
            min_target_idx = target_indices[idx]
            max_target_idx = min_target_idx

            start_idx = max(0, min_target_idx - num_genes)
            end_idx = min(len(loci) - 1, max_target_idx + num_genes)

            start_locus = loci[start_idx]
            end_locus = loci[end_idx]

            start = int(start_locus.location.start)
            end = int(end_locus.location.end)

            subrecord = record[start:end]
            
            # Extract gene name for filename
            gene_name = loci[min_target_idx].qualifiers.get('gene', ['unknown_gene'])[0]
            filename = f"{output_prefix}_{idx+1}_{gene_name}_{record.id}.gbk"
            output_path = os.path.join(output_prefix, filename)
            
            with open(output_path, "w") as output_handle:
                SeqIO.write(subrecord, output_handle, "genbank")
                
            extracted_count += 1
            actual_upstream = min_target_idx - start_idx
            actual_downstream = end_idx - max_target_idx
            
            print(f"Extracted gene region with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id} into {filename}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0    

    return extracted_count

def slice_locustag(record, locus_tags, num_genes, case_insensitive, output_prefix):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    
    try:
        target_indices = [idx for idx, feat in enumerate(loci) if 'locus_tag' in feat.qualifiers and any(locus_tag_matches(feat.qualifiers['locus_tag'][0], tag, case_insensitive) for tag in locus_tags)]
        
        if not target_indices:
            return 0

        extracted_count = 0

        for idx in range(len(target_indices)):
            min_target_idx = target_indices[idx]
            max_target_idx = min_target_idx

            start_idx = max(0, min_target_idx - num_genes)
            end_idx = min(len(loci) - 1, max_target_idx + num_genes)

            start_locus = loci[start_idx]
            end_locus = loci[end_idx]

            start = int(start_locus.location.start)
            end = int(end_locus.location.end)

            subrecord = record[start:end]
            
            # Extract locus tag for filename
            locus_tag = loci[min_target_idx].qualifiers.get('locus_tag', ['unknown_locus_tag'])[0]
            filename = f"{output_prefix}_{idx+1}_{locus_tag}_{record.id}.gbk"
            output_path = os.path.join(output_prefix, filename)
            
            with open(output_path, "w") as output_handle:
                SeqIO.write(subrecord, output_handle, "genbank")
                
            extracted_count += 1
            actual_upstream = min_target_idx - start_idx
            actual_downstream = end_idx - max_target_idx
            
            print(f"Extracted locus tag region with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id} into {filename}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0    

    return extracted_count

def slice_range(record, coords, output_prefix):
    start, end, contig = coords.split(":")
    try:
        start = int(start)
        end = int(end)

        # Locate the specified contig in the record
        for contig in record:
            if contig.name == contig:
                # Slice the specified range within the contig
                subrecord = contig[start:end]
                print(f"Extracted range {start}:{end} from {contig.name}")
                
                # Create a unique filename
                filename = f"{output_prefix}_{record.id}_range_extraction.gbk"
                output_path = os.path.join(output_prefix, filename)
                
                # Write subrecord to file
                with open(output_path, "w") as output_handle:
                    SeqIO.write(subrecord, output_handle, "genbank")
                return subrecord

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def main():
    args = get_args()

    # Check that options are not used simultaneously
    if sum([bool(args.genes), bool(args.range), bool(args.locus_tags)]) > 1:
        sys.stderr.write("Error: -lt, -g and -r options cannot be used simultaneously.\n")
        exit(1)

    records_processed = 0
    records_with_match = 0

    records = SeqIO.parse(args.infile, 'genbank')

    # Use the provided prefix for the output directory
    output_dir = args.outfile
    os.makedirs(output_dir, exist_ok=True)

    for record in records:
        records_processed += 1
        sys.stderr.write(f'Searching {record.id}.\n')
        
        if args.genes:
            extracted_count = slice_gene(record, args.genes.split(':'), args.num_genes, args.case_insensitive, output_dir)
            records_with_match += extracted_count
        
        elif args.locus_tags:
            extracted_count = slice_locustag(record, args.locus_tags.split(':'), args.num_genes, args.case_insensitive, output_dir)
            records_with_match += extracted_count
        
        elif args.range:
            range_parts = args.range.split(':')
            if len(range_parts) != 3:
                sys.stderr.write("Error: When using -r option, provide start:end:contig.\n")
                exit(1)
            start, end, locus = range_parts[0], range_parts[1], range_parts[2]
            if record.name == locus:
                subrecord = slice_range(record, start, end, output_dir)
                if subrecord:
                    # Filename creation based on output prefix and record ID
                    filename = f"{args.outfile}_{record.id}_range_extraction.gb"
                    output_path = os.path.join(output_dir, filename)
                    with open(output_path, "w") as output_handle:
                        SeqIO.write(subrecord, output_handle, "genbank")
                    records_with_match += 1

    # Output summary
    if records_with_match > 0:
        print(f"Found matches in {records_with_match} out of {records_processed} records.")
        print(f"Extracted records are saved in the '{output_dir}' directory.")
    else:
        if args.genes:
            print(f"No matching genes found for {args.genes} in any of the {records_processed} records processed.")
        elif args.locus_tags:
            print(f"No matching locus tags found for {args.locus_tags} in any of the {records_processed} records processed.")
        elif args.range:
            print(f"No matching range found for {args.range} in any of the {records_processed} records processed.")
        else:
            print(f"No matches found in any of the {records_processed} records processed.")

if __name__ == "__main__":
    main()
    
