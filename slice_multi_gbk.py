# Slices/subsections genbank files by either gene names or base pair range. Additionally, can extract target genes ± x number of genes, or ± x bp window. Handles multi-genbank files which have multiple contigs.

# Requirements:
#   - biopython (any)

# Author: Ben Vezina
#   - Scholar: https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en&oi=ao
#   - ORCID: https://orcid.org/0000-0003-4224-2537

# Citation: https://github.com/bananabenana/slice_multi_gbk/
# Adapted from https://gist.github.com/jrjhealey/2df3c65c7a70cbca4862e94620e4d7b2


from Bio import SeqIO
import sys, argparse
import os
import re
from argparse import RawTextHelpFormatter
import fnmatch

def get_args():
    """Parse command line arguments"""
    try:
        parser = argparse.ArgumentParser(
            description='Subset genbanks between 1 gene, 2 genes or base pair ranges.\n'
                        'usage:\n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -o prefix  # Slice between two genes (inclusive of genes) \n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -o prefix  # Slice between two locus tags (inclusive) \n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1:gene2 -n [integer] -o prefix  # Slice between two genes (inclusive) plus [integer] genes up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag1:locus_tag2 -n [integer] -o prefix  # Slice between two locus tags (inclusive) plus [integer] locus tags up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1: -o prefix  # Slice one gene \n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene*: -o prefix  # Slice all genes matching wildcard \n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1: -n [integer] -o prefix  # Slice one gene plus [integer] genes up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag: -n [integer] -o prefix  # Slice one locus tag plus [integer] locus tag up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -g gene1: -w [integer] -o prefix  # Slice one gene plus [integer] basepairs up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -lt locus_tag: -w [integer] -o prefix  # Slice one locus tag plus [integer] basepairs up- and down-stream \n'
                        '  slice_multi_gbk.py -i infile.gbk -r start:end:contig -o prefix  # Slice base pair range and locus/contig. Contig information required to handle multi-genbank files \n',
            formatter_class=RawTextHelpFormatter)

        parser.add_argument('-i', '--infile', action='store', help='Input genbank [.gbk, .gbff] file to slice or subsection. Can be a multi-genbank file with multiple contigs')
        parser.add_argument('-g', '--genes', action='store', help='The two gene names to slice between. Handles wildcards (*). Example: -g hemF:maeB. A single gene can be extracted using one gene name. Example: -g hemF:')
        parser.add_argument('-n', '--num_genes', help='Number of genes to slice before and after the selected genes. Can be used with -g or -lt only. Will take maximum number of genes up/downstream if at end of a contig. Optional. Example: -n 2', type=int, default=0)
        parser.add_argument('-w', '--window', help='Number of base pairs to slice before and after the selected gene/locus tag. Example: -w 10000 for 10kb up/downstream', type=int)
        parser.add_argument('-lt', '--locus_tags', help='The locus tag of the gene to slice if gene name is not available. Handles wildcards (*). Example: -lt DKHHJM_01915:')
        parser.add_argument('-r', '--range', action='store', help='The two base pair coordinates (range) to slice between. Contig name or locus must be provided. Example: 612899:630276:contig_1')
        parser.add_argument('-o', '--outfile', required=True, help='Output directory and filename prefix')
        parser.add_argument('-p', '--protein', action='store_true', help='Produce multifasta files for nucleotide (CDS) and amino acid sequences of all sliced genbanks.')
        parser.add_argument('-c', '--case_insensitive', action='store_true', help='Turn off case sensitivity for gene and locus tag matching.')
        parser.add_argument('-v', '--version', action='version', version='%(prog)s version 3.0.0')

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write('An exception occurred with argument parsing. Check your provided options.')

    return parser.parse_args()

def gene_matches(gene_name, pattern, case_insensitive):
    """Check if a gene name matches the pattern, respecting case sensitivity."""
    print(f"Matching {gene_name} against pattern {pattern}")  # Debugging line
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

    # Iterate over features and collect sequences
    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get('gene', ['unknown_gene'])[0]
            cds_seq = feature.extract(record.seq)
            nucleotides.append(f">{gene_name}\n{cds_seq}\n")
            if 'translation' in feature.qualifiers:
                aa_seq = feature.qualifiers['translation'][0]
                amino_acids.append(f">{gene_name}\n{aa_seq}\n")

    # Ensure output directory exists
    output_dir = os.path.dirname(base_filename)
    os.makedirs(output_dir, exist_ok=True)

    # Write nucleotide multifasta
    if nucleotides:
        fasta_file_fna = f"{base_filename}_nucleotides.fna"
        with open(fasta_file_fna, "w") as fna_file:
            fna_file.writelines(nucleotides)
        print(f"Wrote nucleotide multifasta to {fasta_file_fna}")

    # Write amino acid multifasta
    if amino_acids:
        fasta_file_faa = f"{base_filename}_amino_acids.faa"
        with open(fasta_file_faa, "w") as faa_file:
            faa_file.writelines(amino_acids)
        print(f"Wrote amino acid multifasta to {fasta_file_faa}")


def slice_genomic_data(record, patterns, num_genes, window_size, case_insensitive, output_prefix, produce_multifasta, extraction_counter, is_gene=True):
    # Determine whether we are slicing by gene or locus tag
    if is_gene:
        loci_key = 'gene'
        match_function = gene_matches
        patterns_str = "genes"
    else:
        loci_key = 'locus_tag'
        match_function = locus_tag_matches
        patterns_str = "locus tags"
    
    loci = [feat for feat in record.features if feat.type == "CDS"]
    if not loci:
        print(f"No CDS features found in record {record.id}.")
        return 0

    # Split and clean the patterns
    patterns = [p.strip() for p in patterns]
    patterns = [p for p in patterns if p]

    if not patterns:
        print(f"No {patterns_str} provided.")
        return 0

    start_pattern = patterns[0]
    end_pattern = patterns[1] if len(patterns) >= 2 else None

    # Compile the pattern for wildcard matching
    start_regex = re.compile(start_pattern, re.IGNORECASE if case_insensitive else 0)
    end_regex = re.compile(end_pattern, re.IGNORECASE if case_insensitive else 0) if end_pattern else None

    # Find all genes that match the start pattern
    matching_start_genes = []
    for idx, feat in enumerate(loci):
        if loci_key in feat.qualifiers:
            current_name = feat.qualifiers[loci_key][0]
            if start_regex.match(current_name):
                matching_start_genes.append(idx)

    if not matching_start_genes:
        print(f"No matching {loci_key} found for start pattern '{start_pattern}'")
        return 0

    # Iterate over all matching start genes
    total_extracted = 0
    for start_idx in matching_start_genes:
        # Find the end gene if an end pattern is provided
        end_idx = start_idx  # default to start_idx if no end pattern
        if end_pattern is not None:
            found = False
            for idx in range(start_idx + 1, len(loci)):
                feat = loci[idx]
                if loci_key in feat.qualifiers:
                    current_name = feat.qualifiers[loci_key][0]
                    if end_regex.match(current_name):
                        end_idx = idx
                        found = True
                        break
            if not found:
                print(f"No matching {loci_key} found for end pattern '{end_pattern}' after start pattern '{start_pattern}'")
                continue  # skip to the next start gene

        # Apply num_genes if provided
        if num_genes > 0:
            new_start_idx = max(0, start_idx - num_genes)
            new_end_idx = min(len(loci) - 1, end_idx + num_genes)
            start_idx = new_start_idx
            end_idx = new_end_idx

        start_locus = loci[start_idx]
        end_locus = loci[end_idx]

        start_pos = int(start_locus.location.start)
        end_pos = int(end_locus.location.end)

        # Apply window if provided
        if window_size is not None:
            desired_start = start_pos - window_size
            desired_end = end_pos + window_size
            start_pos = max(0, desired_start)
            end_pos = min(len(record.seq), desired_end)

        # Ensure the output directory exists
        os.makedirs(output_prefix, exist_ok=True)

        window_suffix = f"_{window_size}bp_window" if window_size is not None else ""
        num_genes_suffix = f"_{num_genes}genes" if num_genes > 0 else ""
        # Include start_idx and end_idx in the filename to make it unique
        base_filename = os.path.join(output_prefix, f"{output_prefix}_{extraction_counter}_{start_idx}_{end_idx}_{start_pattern}_to_{end_pattern if end_pattern else ''}{num_genes_suffix}{window_suffix}_{record.id}")

        # Write GBK file
        filename_gbk = f"{base_filename}.gbk"
        with open(filename_gbk, "w") as gbk_file:
            subrecord = record[start_pos:end_pos]
            SeqIO.write(subrecord, gbk_file, "genbank")

        print(f"Extracted {patterns_str} matching patterns {patterns} from {start_pos} to {end_pos} into {filename_gbk}")

        # Write multifasta files if requested
        if produce_multifasta:
            write_multifasta(subrecord, base_filename, produce_multifasta)

        total_extracted += 1

    return total_extracted


def slice_range(record, start, end, contig, output_prefix, produce_multifasta):
    try:
        if contig == record.id:
            subrecord = record[start:end]
            
            filename = f"{output_prefix}_{start}_{end}_{contig}.gbk"
            output_path = os.path.join(output_prefix, filename)

            # Ensure the output directory exists
            os.makedirs(output_prefix, exist_ok=True)

            with open(output_path, "w") as output_handle:
                SeqIO.write(subrecord, output_handle, "genbank")
                
            print(f"Extracted range from {start} to {end} on contig {contig} into {filename}")
            
            # Write multifasta files if requested
            if produce_multifasta:
                # Create a base filename for multifasta files
                base_filename = os.path.join(output_prefix, f"{output_prefix}_{start}_{end}_{contig}")
                write_multifasta(subrecord, base_filename, produce_multifasta)

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return 0

    return 1
def main():
    args = get_args()
    input_file = args.infile
    output_prefix = args.outfile
    produce_multifasta = args.protein
    case_insensitive = args.case_insensitive
    window_size = args.window

    extraction_counter = 1

    try:
        os.makedirs(output_prefix, exist_ok=True)

        with open(input_file, "r") as file_handle:
            records = SeqIO.parse(file_handle, "genbank")
            
            total_extracted = 0

            for record in records:
                if args.genes:
                    gene_list = args.genes.split(':')
                    extracted = slice_genomic_data(record, gene_list, args.num_genes, window_size, case_insensitive, output_prefix, produce_multifasta, extraction_counter, is_gene=True)
                    if extracted:
                        total_extracted += extracted
                        extraction_counter += 1

                if args.locus_tags:
                    locus_tag_list = args.locus_tags.split(':')
                    extracted = slice_genomic_data(record, locus_tag_list, args.num_genes, window_size, case_insensitive, output_prefix, produce_multifasta, extraction_counter, is_gene=False)
                    if extracted:
                        total_extracted += extracted
                        extraction_counter += 1

                if args.range:
                    range_components = args.range.split(':')
                    if len(range_components) != 3:
                        print("Invalid range format. Expected start:end:contig")
                        continue
                    range_start, range_end, contig = range_components
                    extracted = slice_range(record, int(range_start), int(range_end), contig, output_prefix, produce_multifasta)
                    if extracted:
                        total_extracted += extracted

            print(f"Total records extracted: {total_extracted}")
    
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
