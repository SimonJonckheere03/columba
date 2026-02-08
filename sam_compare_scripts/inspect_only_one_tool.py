"""
This script inspects SAM files named with the prefix "only_" (excluding those containing  "unknown") in a given directory.

Each input SAM file lacks a header, which is provided separately and used to parse the alignments correctly.

For each unique read ID, the script evaluates its alignments and categorizes them based on:

- Presence of soft (S) or hard (H) clipping in the CIGAR string,
- Alignment error rate exceeding 8% of read length (using the NM tag),
- Presence of non-ACGT characters in the reference sequence,
- Read length shorter than 100 bases,
- Or marked unknown if none of these conditions apply.

It outputs counts of unique reads and the number in each category, and writes read IDs with unknown classification to separate files for further review.

It is intented to be run as a part of the compare_sam.py script and generally does not need to be run separately.

Usage:
    python inspect_missed_Columba.py <input_directory> <sam_file_with_header>
"""
import sys
import os
from Bio import SeqIO
import pysam

from utils import  s_or_h_in_cigar, contains_non_acgt
from extract_sub_from_fasta import load_sequences



def read_sam_with_external_header(file, header_obj):
    """
    Open a file with an external header.
    """
    lines = []
    with open(file, 'r') as f:
        for line in f:
            aln = pysam.AlignedSegment.fromstring(line, header_obj)
            lines.append(aln)
    return lines
                               
def read_length(cigar_tuples):
    """
    Calculate the length of the read based on the CIGAR string.
    The read length is defined as the sum of lengths of M, I, S, and H operations.
    """
    return sum(count for op, count in cigar_tuples if op in (0, 1, 4, 5))  # M, I, S or H

def error_rate(score, cigar_tuples):
    """
    Calculate the error rate based on the CIGAR string.
    The error rate is defined as the number of mismatches divided by the length of the aligned sequence.
    """
    length = read_length(cigar_tuples)
    if length == 0:
        return 0.0
    return score / length

def error_rate_record(alignmentRecord):
    """
    Calculate the error rate based on the CIGAR string.
    The error rate is defined as the number of mismatches divided by the length of the aligned sequence.
    """
    score = alignmentRecord.get_tag("NM")  # Assuming NM tag is used for mismatch count
    if score is None:
        return 0.0
    return error_rate(score, alignmentRecord.cigar)


def process_read_id(read_id, numbers_dict, alignments, sequences, unknown_reasons, threshold=0.08):
    """
    Process a read ID and update the numbers dictionary based on the alignments.
    """
    numbers_dict["unique_read_ids"] += 1

    if all(s_or_h_in_cigar(aln.cigar) for aln in alignments):
        numbers_dict["SorH"] += 1
    elif all(error_rate_record(aln) > threshold for aln in alignments):
        numbers_dict["error_rate_above_threshold"] += 1
    else:
        aln_format = [(aln.reference_name, aln.reference_start, aln.cigartuples) for aln in alignments]
        if contains_non_acgt(aln_format, sequences):
            numbers_dict["non_acgt_characters"] += 1
        elif read_length(alignments[0].cigar) < 100:  
            numbers_dict["short read"] += 1
        else:
            numbers_dict["unknown_reason"] += 1
            unknown_reasons.append((read_id, alignments))



    return numbers_dict, unknown_reasons

def inspect_only(only_file, header_obj, sequences, threshold=0.08):
    """
    Inspect the only file for unique read IDs and their alignments.
    """
    sam_file = read_sam_with_external_header(only_file, header_obj)
    
    current_read_id = None
    alignments_current_read_id = []

    numbers_dict = {
        "unique_read_ids": 0,
        "SorH": 0,
        "error_rate_above_threshold": 0,
        "unknown_reason": 0,
        "short read": 0,
        "non_acgt_characters": 0,
        "count": 0}
    
    unknown_reasons = []
    
    for alignment in sam_file:
        read_id = alignment.query_name
        if read_id == current_read_id:
            alignments_current_read_id.append(alignment)
            continue
        
        # Process the previous read ID
        if current_read_id is not None:
            numbers_dict, unknown_reasons = process_read_id(current_read_id, numbers_dict, alignments_current_read_id, sequences, unknown_reasons, threshold)

        # Reset for the new read ID
        current_read_id = read_id
        alignments_current_read_id = [alignment]

    # Process the last read ID
    if current_read_id is not None:
        numbers_dict, unknown_reasons = process_read_id(current_read_id, numbers_dict, alignments_current_read_id, sequences, unknown_reasons,  threshold)

    numbers_dict["count"] = numbers_dict["unique_read_ids"]
    return numbers_dict, unknown_reasons


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python inspect_missed_Columba.py <input_directory> <sam_file_with_header>")
        sys.exit(1)

    input_directory = sys.argv[1]
    sam_file_with_header = sys.argv[2]
    header_obj = pysam.AlignmentFile(sam_file_with_header, "r").header

    # list files starting with "only_" but that do not contain "Columba" in their name
    only_files = [f for f in os.listdir(input_directory) if f.startswith("only_")  and "unknown" not in f]

    threshold = 0.08  # 8% error rate threshold

    fasta_file = "/data/lrenders/newIndex/hs.grch38.fa"
    sequences = load_sequences(fasta_file)  # Load sequences from the FASTA file

    # Open the SAM file with the external header
    for f in only_files:
        # if file ends in ".txt", skip it
        if f.endswith(".txt"):
            continue

        numbers_dict, unknown_reasons = inspect_only(os.path.join(input_directory, f), header_obj, sequences, threshold)
        unknown_file = open(os.path.join(input_directory, f"{f}_unknown_reasons.txt"), "w")


        # write unknown reasons to file
        for read_id, alignments in unknown_reasons:
            for aln in alignments:
                unknown_file.write(f"{aln.query_name}\n")
        unknown_file.close()
 
        # Print the results
        print(f"Results for {f}:")
        print("Number of unique read IDs:", numbers_dict["unique_read_ids"])
        print("Number of reads with S or H in CIGAR string:", numbers_dict["SorH"])
        print("Number of reads with error rate above 8%:", numbers_dict["error_rate_above_threshold"])
        print("Number of reads with unknown reasons:", numbers_dict["unknown_reason"])
        print(f"\tOf which {numbers_dict['non-acgt_handling_unknown_reason']} have non-ACGT characters in the reference sequence")       
                    



