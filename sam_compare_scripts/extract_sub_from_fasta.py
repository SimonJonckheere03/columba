from Bio import SeqIO
import sys
import gzip

"""
Extract a subsequence from a FASTA file based on sequence name, start position, and length.
It assumes the start position is 1-based indexing, (like in SAM format).
Usage:
    python extract_sub_from_fasta.py <seq_name> <start_position> <length> <fasta_file>
Additonally, this file provides functions to extract multiple subsequences from a dictionary of sequences
and to list all sequences in a FASTA file (for repeated use).

"""

def extract_subsequence(fasta_file, seq_name, start_position, length):
    # Check if the file is compressed
    if fasta_file.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"  # 't' stands for text mode (crucial for SeqIO)
    else:
        open_func = open
        mode = "r"

    # Open and parse the FASTA file using the selected method
    with open_func(fasta_file, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == seq_name:
                # Since start_position is 1-based, convert to 0-based indexing
                start = start_position - 1
                end = start + length
                # Extract the substring
                return str(record.seq[start:end])
    
    return None  # Return None if sequence not found


def load_sequences(fasta_file):
    # Determine how to open the file based on extension
    if fasta_file.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"  # Read Text mode is required for gzip
    else:
        open_func = open
        mode = "r"

    # Open the file and parse
    with open_func(fasta_file, mode) as handle:
        return {record.id: record.seq for record in SeqIO.parse(handle, "fasta")}

def extract_multiple_subsequences(sequences, queries):
    """ Extract multiple subsequences from a dictionary of sequences based on queries.
    Args:
        sequences (dict): A dictionary where keys are sequence names and values are Seq objects.
        queries (list): A list of tuples, each containing (seq_name, start_position, length).
    Returns:
        list: A list of extracted subsequences as strings.
    """
    results = []

    # sort queries by seq_name and start_position
    queries.sort(key=lambda x: (x[0], x[1]))
    for seq_name, start_position, length in queries:
        seq = sequences.get(seq_name)
        if seq:
            start = start_position - 1
            results.append(str(seq[start:start + length]))
        else:
            results.append(None)
    return results

if __name__ == "__main__":

    # name, pos and length are arguments
    args = sys.argv[1:]
  
    if   len(args) != 4:
        print("Usage: extract_sub_from_fasta.py seq_name start_position length fasta_file")
        sys.exit(1)

    seq_name = args[0]
    start_position = int(args[1])
    length = int(args[2])
    fasta_file = args[3]


    subsequence = extract_subsequence(fasta_file, seq_name, start_position, length)
    if subsequence:
        print(subsequence)
    else:
        print(f"Sequence '{seq_name}' not found in '{fasta_file}'")
