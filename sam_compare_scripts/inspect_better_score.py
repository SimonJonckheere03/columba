"""
This script analyzes files named score_better_tool.sam to investigate cases where one alignment tool
produces a better score than another. 

Input files have the following structure per read:
1. Read ID line
2. Info line for tool 1: number of alignments
3. Corresponding alignment lines for tool 1
4. Info line for tool 2: number of alignments
5. Corresponding alignment lines for tool 2
6. Blank line separating entries

The script parses these files, compares alignments between the two tools, and categorizes differences 
based on criteria such as presence of soft/hard clipping in CIGAR strings, positional similarity, and
presence of non-ACGT characters in the aligned sequences.

It is intended as a helper for compare_sam.py and assumes that compare_sam.py's initial processing
step has already been performed.
It is ran as a part of compare_sam.py and generally does not need to be run separately.

Usage:
    python inspect_better_score.py <input_directory> <sam_file_with_header> <edit_distance>
"""


import pysam
import os
import sys



from extract_sub_from_fasta import load_sequences  
from utils import contains_non_acgt, make_dict, match_other_as_dict, s_or_h_in_cigar 




def parse_alignment_block(file, n, header):
    alignments = []
    for _ in range(n):
        line = file.readline().strip()
        aln = pysam.AlignedSegment.fromstring(line, header)
        alignments.append((aln.reference_name, aln.reference_start, aln.cigartuples))
    return alignments


    # Compare

    
def match(pos, lst, tolerance):
    """
    Check if the position matches any position in the list within a given tolerance.
    """
    return any(r == pos[0] and abs(rp - pos[1]) <= tolerance for r, rp, _ in lst)




def parse_file(path, header, tolerance, sequences):
    
    # loop over lines in the file
    current_read_id = None

    # read the first line
    f = open(path, "r")
    next_read_id = f.readline().strip()

    numbers_dict= {"SorH": 0, "relative_same_positions": 0, "non_acgt_handling_rel_same_pos": 0, "widely_different_alignment": 0,
                   "non_acgt_handling_widely_different": 0, "unknown_reasons": 0, "count": 0}

    relative_same_positions = []
    unknown_reasons = []

    processed = 0

    while next_read_id:
        # Read the next read ID
        current_read_id = next_read_id
        numbers_dict["count"] += 1

        # Read the info line for tool 1
        info_1 = f.readline().strip()
        n_1 = int(info_1.split(":")[1].split()[0])
        alignments_1 = parse_alignment_block(f, n_1, header)

        # Read the info line for tool 2
        info_2 = f.readline().strip()
        n_2 = int(info_2.split(":")[1].split()[0])
        alignments_2 = parse_alignment_block(f, n_2, header)

        # Skip the blank line
        f.readline()

        # find the reason for better score
        # option 1: S or H in CIGAR string

        isSorH = any(s_or_h_in_cigar(a[2]) for a in alignments_1) or any(s_or_h_in_cigar(a[2]) for a in alignments_2)
      
        numbers_dict["SorH"] += isSorH
        if not isSorH:
            # check if the alignments are on the same reference and similar positions

            alignments_1_dict = make_dict(alignments_1)
            alignments_2_dict = make_dict(alignments_2)
            any_simlar = any(match_other_as_dict(a, alignments_2_dict, tolerance) for a in alignments_1) and \
            any(match_other_as_dict(a, alignments_1_dict, tolerance) for a in alignments_2)

            if any_simlar:
                numbers_dict["relative_same_positions"] += 1
                relative_same_positions.append((current_read_id, alignments_1, alignments_2))
            else:
                numbers_dict["widely_different_alignment"] += 1
                unknown_reasons.append((current_read_id, alignments_1, alignments_2))
        # Read the next read ID
        next_read_id = f.readline().strip()

    relative_same_pos_no_acgt= []
    unknown_reasons_no_acgt = []
    if numbers_dict["relative_same_positions"] > 0 or numbers_dict["widely_different_alignment"] > 0:
     
        # check if the aligned sequences of one of the tools contains a non-ACGT character
        for read_id, a1, a2 in relative_same_positions:
            if contains_non_acgt(a1, sequences) or contains_non_acgt(a2, sequences):
                numbers_dict["non_acgt_handling_rel_same_pos"] += 1
            else:
                relative_same_pos_no_acgt.append((read_id, a1, a2))
                

        for read_id, a1, a2 in unknown_reasons:
            if contains_non_acgt(a1, sequences) or contains_non_acgt(a2, sequences):
                numbers_dict["non_acgt_handling_widely_different"] += 1
            else:
                unknown_reasons_no_acgt.append((read_id, a1, a2))
                


    return unknown_reasons_no_acgt, relative_same_pos_no_acgt, numbers_dict

    

# Usage
if __name__ == "__main__":
    
    # get input directory from user
    if len(sys.argv) != 4:
        print("Usage: python inspect_better_score.py <input_directory> <sam_file_with_header> <edit distance>")
        sys.exit(1)

    input_directory = sys.argv[1]
    sam_file_with_header = sys.argv[2]

    ref_sam = pysam.AlignmentFile(sam_file_with_header, "r")  # or "rb" for BAM
    header_obj = ref_sam.header  # this is already a pysam.AlignmentHeader
    edit_distance = int(sys.argv[3]) 

    tolerance = 2* edit_distance + 1  # tolerance for matching positions

    fasta_file ="/data/lrenders/newIndex/hs.grch38.fa" # todo make this a parameter

    sequences = load_sequences(fasta_file)  # Load sequences from the FASTA file
    #do this for all files in directory starting with "score_better_"

    file_score_better_ids= open(os.path.join(input_directory, "score_better_ids.txt"), "w")

    for filename in os.listdir(input_directory):

        if filename.startswith("score_better_"):
            file_score_better_ids.write(filename + "\n")
            results = parse_file(os.path.join(input_directory, filename), header_obj, tolerance, sequences)
            unknown_reasons, relative_same_positions, numbers_dict = results
            print(f"Results for {filename}:")
            print("Number of reads with S or H in CIGAR string:", numbers_dict["SorH"])
            print("Number of reads with relative same positions:", numbers_dict["relative_same_positions"])
            print("\tOf which had non-ACGT characters in the aligned sequences of one of the tools:", numbers_dict["non_acgt_handling_rel_same_pos"])
            print("Number of reads with Widely different alignments:", numbers_dict["widely_different_alignment"])
            print("\tOf which had non-ACGT characters in the aligned sequences of one of the tools:", numbers_dict["non_acgt_handling_widely_different"])
      
    
            print("Relative same positions but no non-ACGT character in reference:")
            file_score_better_ids.write("Relative same positions but no non-ACGT character in reference:\n")
            for read_id, a1, a2 in  relative_same_positions:
                print(f"Read ID: {read_id}")
                print("Tool1:")
                for seq_id, pos, cigar in a1:
                    print(f"\t{seq_id} {pos} {cigar}")
                print("Tool2:")
                for seq_id, pos, cigar in a2:
                    print(f"\t{seq_id} {pos} {cigar}")
                file_score_better_ids.write(f"{read_id}\n")

  
            print("Unknown reasons:")
            file_score_better_ids.write("Unknown reasons (no non-ACGT):\n")
            for read_id, a1, a2 in unknown_reasons:
                print(f"Read ID: {read_id}")
                print("Tool1:")
                for seq_id, pos, cigar in a1:
                    print(f"\t{seq_id} {pos} {cigar}")
                print("Tool2:")
                for seq_id, pos, cigar in a2:
                    print(f"\t{seq_id} {pos} {cigar}")
                file_score_better_ids.write(f"{read_id}\n")
            input("Press Enter to continue...")
   
