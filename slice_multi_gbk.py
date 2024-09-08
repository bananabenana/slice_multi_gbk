#!/usr/bin/env python3
# Slices/subsections genbank files by either gene names or base pair range. This handles multi-genbank files which have multiple contigs

# Usage:
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output # This will extract the region between two genes including these genes
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o gene_sliced_output -p # This will extract the region between two genes including these genes, also outputting multi-fasta protein and CDS files too
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o locus_tag_sliced_output # This will extract the region between two locus tags if gene names not available
#   slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n 2 -o gene_sliced_output # This will extract 2 genes up- and down-stream of your specified two genes
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n 2 -o locus_tag_sliced_output
#   slice_multi_gbk.py -i infile.gbk -g gene1: -o gene_extract_output # This will extract a single gene
#   slice_multi_gbk.py -i infile.gbk -g gene*: -o wildcard_gene_extract_output # This will extract any genes matching this wildcard
#   slice_multi_gbk.py -i infile.gbk -g gene1: -n 3 -o gene_extract_output # This will extract a single gene and 3 genes up- and down-stream of the gene
#   slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n 3 -o gene_extract_output # This will extract a single locus tag and 3 genes up- and down-stream of this
#   slice_multi_gbk.py -i infile.gbk -r start:end:contig -o range_sliced_output # This will slice via base pair range and locus/contig. This is required to handle multi-genbank files

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

#!/usr/bin/env python3
# Slices/subsections genbank files by either gene names or base pair range. This handles multi-genbank files which have multiple contigs

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
        parser.add_argument('-p', '--protein', action='store_true', help='Produce multifasta files for nucleotide (CDS) and amino acid sequences of all sliced genbanks.')
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

def write_multifasta(record, base_filename, produce_multifasta):
    """Write multifasta files for nucleotide (CDS) and amino acid sequences."""
    if not produce_multifasta:
        return

    nucleotides = []
    amino_acids = []

    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get('gene', ['unknown_gene'])[0]
            cds_seq = feature.extract(record.seq)
            nucleotides.append(f">{gene_name}\n{cds_seq}\n")
            if 'translation' in feature.qualifiers:
                aa_seq = feature.qualifiers['translation'][0]
                amino_acids.append(f">{gene_name}\n{aa_seq}\n")

    # Write nucleotide multifasta
    if nucleotides:
        fasta_file_fna = f"{base_filename}.fna"
        with open(fasta_file_fna, "w") as fna_file:
            fna_file.writelines(nucleotides)
        print(f"Wrote nucleotide multifasta to {fasta_file_fna}")

    # Write amino acid multifasta
    if amino_acids:
        fasta_file_faa = f"{base_filename}.faa"
        with open(fasta_file_faa, "w") as faa_file:
            faa_file.writelines(amino_acids)
        print(f"Wrote amino acid multifasta to {fasta_file_faa}")

def process_records(input_file, genes, num_genes, output_prefix, produce_multifasta):
    extraction_counter = 1  # Initialize extraction counter

    # Parse the input file
    for record in SeqIO.parse(input_file, "genbank"):
        # Call the function and pass the extraction_counter
        success = slice_gene(record, genes, num_genes, output_prefix, produce_multifasta, extraction_counter)
        
        # Increment the extraction_counter
        if success:
            extraction_counter += 1

def slice_gene(record, genes_patterns, num_genes, output_prefix, produce_multifasta, extraction_counter, case_insensitive):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    
    try:
        # Find indices of the target genes
        target_indices = [
            idx for idx, feat in enumerate(loci)
            if 'gene' in feat.qualifiers and any(gene_matches(feat.qualifiers['gene'][0], pattern, case_insensitive) for pattern in genes_patterns)
        ]

        if not target_indices:
            print(f"No matching genes found for patterns: {genes_patterns}")
            return 0

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

        # Ensure the output directory exists
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)

        # Generate base filename with extraction counter
        base_filename = os.path.join(output_prefix, f"{output_prefix}_{extraction_counter}_{genes_patterns[0]}_{record.id}")

        # Write GBK file
        filename_gbk = f"{base_filename}.gbk"
        with open(filename_gbk, "w") as gbk_file:
            SeqIO.write(subrecord, gbk_file, "genbank")

        print(f"Extracted genes matching patterns {genes_patterns} with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id} into {filename_gbk}")

        # Write multifasta files if requested
        if produce_multifasta:
            write_multifasta(subrecord, base_filename, produce_multifasta)

        return 1

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0

def slice_locustag(record, locus_tags_patterns, num_genes, case_insensitive, output_prefix, produce_multifasta, extraction_counter):
    loci = [feat for feat in record.features if feat.type == "CDS"]
    
    try:
        # Find indices of the target locus tags
        target_indices = [
            idx for idx, feat in enumerate(loci)
            if 'locus_tag' in feat.qualifiers and any(locus_tag_matches(feat.qualifiers['locus_tag'][0], pattern, case_insensitive) for pattern in locus_tags_patterns)
        ]

        if not target_indices:
            print(f"No matching locus tags found for patterns: {locus_tags_patterns}")
            return 0

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

        # Extract the subrecord for the entire region
        subrecord = record[start:end]

        actual_upstream = min_target_idx - start_idx
        actual_downstream = end_idx - max_target_idx

        # Ensure the output directory exists
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)

        # Generate base filename with extraction counter
        base_filename = os.path.join(output_prefix, f"{output_prefix}_{extraction_counter}_locus_tags_{record.id}")

        # Write GBK file
        filename_gbk = f"{base_filename}.gbk"
        with open(filename_gbk, "w") as output_handle:
            SeqIO.write(subrecord, output_handle, "genbank")

        print(f"Extracted locus tag regions with {actual_upstream} gene(s) upstream and {actual_downstream} gene(s) downstream from {record.id} into {filename_gbk}")

        # Write multifasta files if requested
        if produce_multifasta:
            write_multifasta(subrecord, base_filename, produce_multifasta)

        return 1

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0

def slice_range(record, start, end, contig, output_prefix, produce_multifasta):
    try:
        if contig == record.id:
            subrecord = record[start:end]
            
            filename = f"{output_prefix}_{contig}_{start}_{end}.gbk"
            output_path = os.path.join(output_prefix, filename)
            
            with open(output_path, "w") as output_handle:
                SeqIO.write(subrecord, output_handle, "genbank")
                
            print(f"Extracted range from {start} to {end} on contig {contig} into {filename}")
            
            if produce_multifasta:
                write_multifasta(subrecord, output_prefix, produce_multifasta)
                    
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0

    return 1

def main():
    args = get_args()
    input_file = args.infile
    output_prefix = args.outfile
    produce_multifasta = args.produce_multifasta
    case_insensitive = args.case_insensitive  # Retrieve the case insensitive flag

    try:
        # Make output directory if it does not exist
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)

        # Initialize extraction counter
        extraction_counter = 1

        with open(input_file, "r") as file_handle:
            records = SeqIO.parse(file_handle, "genbank")
            
            total_extracted = 0

            for record in records:
                if args.genes:
                    gene_list = args.genes.split(':')
                    extracted = slice_gene(record, gene_list, args.num_genes, output_prefix, produce_multifasta, extraction_counter, case_insensitive)
                    total_extracted += extracted
                    extraction_counter += 1  # Increment extraction counter after each slice

                if args.locus_tags:
                    locus_tag_list = args.locus_tags.split(':')
                    extracted = slice_locustag(record, locus_tag_list, args.num_genes, case_insensitive, output_prefix, produce_multifasta, extraction_counter)
                    total_extracted += extracted
                    extraction_counter += 1  # Increment extraction counter after each slice

                if args.range:
                    range_start, range_end, contig = args.range.split(':')
                    extracted = slice_range(record, int(range_start), int(range_end), contig, output_prefix, produce_multifasta)
                    total_extracted += extracted

            print(f"Total records extracted: {total_extracted}")
    
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
