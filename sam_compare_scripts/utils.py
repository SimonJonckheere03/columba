# Utility functions for working with alignments (typically from SAM files)
# Includes:
# - Detection of non-ACGT characters in aligned reference regions
# - CIGAR string interpretation for alignment and read lengths
# - Clipping checks (S/H)
# - Efficient positional matching with tolerance using binary search
# - Conversion of alignment tuples into a reference-position indexed dictionary
#
# Assumes alignment tuples follow the format:
# (reference_name, position, cigar_tuples, [edit_distance])
#
# The 'contains_non_acgt' function depends on extract_sub_from_fasta.extract_multiple_subsequences
# and a dictionary of reference sequences.


from extract_sub_from_fasta import extract_multiple_subsequences 
from bisect import bisect_left
from collections import defaultdict

def contains_non_acgt(alignments, sequences):
    """
    Check if the aligned sequences contains a non-ACGT character.
    Returns: True if any of the aligned sequences contains a non-ACGT character, False otherwise.
    alignments: List of tuples (seq_name, start_position, cigarstring)
    sequences: Dictionary of sequences where keys are sequence names and values are the sequences.
    """
    queries = [(seq_name, start_position + 1, alignment_length(cigarstring) + 1 ) for seq_name, start_position, cigarstring in alignments]
    subsequences = extract_multiple_subsequences(sequences, queries)
 
    return any(any(c not in 'ACGT' for c in subsequence) for subsequence in subsequences)


def match_other_as_dict(pos, other, tolerance):
    # other is a dict with seq_id as key and and a list of tuples (pos, cigar_tuples, edit_distance) sorted on pos
    if pos[0] not in other:
        return False
    # since other[pos[0]] is a sorted list, we can use binary search to find the position
    
    index = bisect_left(other[pos[0]], (pos[1],))
    if index < len(other[pos[0]]) and abs(other[pos[0]][index][0] - pos[1]) <= tolerance:
        return True
    if index > 0 and abs(other[pos[0]][index - 1][0] - pos[1]) <= tolerance:
        return True
    if index + 1 < len(other[pos[0]]) and abs(other[pos[0]][index + 1][0] - pos[1]) <= tolerance:
        return True
    return False

def alignment_length(cigartuples):
    """
    Calculate the length of the aligned sequence based on the CIGAR string.
    """
    length = 0
    for op, count in cigartuples:
        if op in (0, 2):  # M, D
            length += count
    return length


def read_length(cigartuples):
    """
    Calculate the length of the read based on the CIGAR string.
    """
    length = 0
    for op, count in cigartuples:
        if op in (0, 1, 4, 5):  # M, I, S, H
            length += count
    return length

def s_or_h_in_cigar(tuples):
    """
    Check if the CIGAR string of the alignment contains 'S' or 'H'.
    These are operations 4 (soft clipping) and 5 (hard clipping) in CIGAR notation.
    """
    return  any(op in (4, 5) for op, _ in tuples)

def make_dict(alignments):
    """
    Make a dictionary from the alignments, where the key is the reference name and the value is a list of tuples (position, cigar_tuples, edit_distance)
    """

    alignments_dict = defaultdict(list)
    for aln in alignments:
        alignments_dict[aln[0]].append(tuple(aln[1:]))
    # sort all key, values in alignments_dict by position
    for key in alignments_dict:    
        alignments_dict[key].sort(key=lambda x: x[0])
 
    return alignments_dict