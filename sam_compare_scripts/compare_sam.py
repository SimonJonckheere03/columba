"""
This file compares two SAM files and finds the differences in alignments.
It first cleans the SAM files by updating CIGAR strings (if none was present), 
removing duplicates, splitting up addititonal alignments in the XA:tag and keeping 
only the best NM alignments.
It checks for:
1. Reads that are only present in one of the files.
2. Reads that have different NM scores (edit distances).
3. Reads that have the same NM score but different number of alignments.
4. Reads that have the same NM score and number of alignments but different reference positions.
The results are written to files in a directory named "compare_<directory_name>".
A detailed analysis of the differences is written to standard output.

Usage:
      python script.py <file1> <file2> <fasta_file> [edit_distance] [threshold] [number_reads]

    Arguments:
      <file1>        : Path to the first SAM file
      <file2>        : Path to the second SAM file
      <fasta_file>   : Path to the reference FASTA file
      [edit_distance]: The maximla edit distance for the lossless tools. Optional integer, default=12, must be non-negative
      [threshold]    : The maximal error rate for the lossless tools. Optional float, default=0.08
      [number_reads] : Optional integer, default=1,000,000
"""
import pysam
import os
import time
import random

import inspect_more_alignments
import inspect_better_score
from inspect_only_one_tool import inspect_only
from extract_sub_from_fasta import load_sequences


# Get the current working directory
current_directory = os.getcwd()


def update_cigar(sam_file, output_file):
    '''
    If no CIGAR provided add XM with X being the length of the query sequence as CIGAR string
    Also add a tag to the read indicating that the CIGAR was updated
    '''
    print("Updating CIGAR strings in file: ", sam_file)
    # Suppress warnings
    original_verb = pysam.get_verbosity()

    i = 0

    pysam.set_verbosity(0)
    # Open the input SAM file
    with pysam.AlignmentFile(sam_file, "r") as sam_in:
        # Create an output SAM file
        with pysam.AlignmentFile(output_file, "w", header=sam_in.header) as sam_out:
            updated_reads = []
            for read in sam_in:
                if read.is_unmapped and read.has_tag('NM'):

                    read.is_unmapped = False  # mark the read as mapped

                    # create a dummy cigar string
                    read.cigar = [(0, read.query_length)]  # 0 = match in CIGAR
                    # indicate that the CIGAR was updated
                    read.set_tag("YC", True)
                updated_reads.append(read)

                if read.has_tag('NM') and read.is_unmapped:
                    print("Updating CIGAR string failed for read: ",
                          read.query_name)
                    print(read.cigarstring)
                    exit(1)

                if len(updated_reads) == 500000:
                    for read in updated_reads:
                        sam_out.write(read)
                    updated_reads = []
                    i += 1
                    print("Processed ", i*500000, " alignments")

            for read in updated_reads:
                sam_out.write(read)

    pysam.set_verbosity(original_verb)


def sort_by_query(input, output):
    pysam.sort("-n", "-o", output, input)


def remove_duplicates(input_file, output_file):
    '''
    Removes duplicate lines from a file. Two lines are considered duplicates if:
    1. The entire line matches, or
    2. All columns except the second match, and the second column differs by 256.
    '''
    def is_duplicate256(existing_line, new_line):
        # Split the lines into columns
        existing_parts = existing_line.split()
        new_parts = new_line.split()

        if len(existing_parts) < 6:
            print("Invalid line: ", existing_line)
        if len(new_parts) < 6:
            print("Invalid line: ", new_line)

        # Check if all columns except the second match
        if existing_parts[0] == new_parts[0] and existing_parts[2:6] == new_parts[2:6]:
            # Check if the second column differs by exactly 256
            return abs(int(existing_parts[1]) - int(new_parts[1])) == 256

        return False

    lines_to_write = []

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        previous_first_column = None
        lines_removed = 0
        count = 0

        current_set = set()
        current_count = 0

        for current_line in infile:
            count += 1
            current_line = current_line.strip()  # Remove trailing whitespace

            if current_line[0] == '@':
                lines_to_write.append(current_line)
                continue

            parts = current_line.split()
            first_column_value = parts[0]

            if first_column_value == previous_first_column:
                # add to the set
                current_set.add(current_line)
                current_count += 1

            else:
                # Write the unique lines of the previous group to the output file
                if previous_first_column is not None:
                    unique_lines = []
                    zeroline = None
                    for line in current_set:
                        # check if second column is 256 or 272
                        if (int(line.split()[1]) != 256 and int(line.split()[1]) != 272):
                            zeroline = line
                        else:
                            unique_lines.append(line)
                        # check if zero line is present in unique lines
                    if zeroline is not None:
                        if not any(is_duplicate256(existing_line, zeroline) for existing_line in unique_lines):
                            unique_lines.append(zeroline)

                    lines_to_write.extend(unique_lines)
                    lines_removed += current_count - len(unique_lines)
                    ID = unique_lines[0].split()[0]
                    if current_count != len(unique_lines):
                        print(
                            f"Removed {current_count - len(unique_lines)} duplicate lines for ID {ID}")
                        print(len(current_set), len(
                            unique_lines), current_count)

                # Start a new group
                current_count = 1
                current_set = set()
                current_set.add(current_line)

                # check if lines_to_write is too large
                if len(lines_to_write) > 1000000:
                    outfile.write('\n'.join(lines_to_write) + "\n")
                    lines_to_write = []

            # Update the previous first column value
            previous_first_column = first_column_value

        # Write the last group of unique lines to the output file
        if current_set:
            unique_lines = []
            for line in current_set:
                # find if there is any line in set for which is_duplicate256 is true
                if not any(is_duplicate256(existing_line, line) for existing_line in unique_lines):
                    unique_lines.append(line)

            lines_to_write.extend(unique_lines)
            lines_removed += current_count - len(unique_lines)
        outfile.write('\n'.join(lines_to_write) + "\n")

    print(
        f'Removed {lines_removed} duplicate lines from {input_file} and wrote to {output_file}')


def remove_unmapped(input_file, output_file):
    """Removes unmapped reads from a SAM/BAM file."""
    write_reads = []
    with pysam.AlignmentFile(input_file, "r") as input_sam:
        with pysam.AlignmentFile(output_file, "w", header=input_sam.header) as output_sam:
            for read in input_sam:
                if not read.is_unmapped:  # Check if the read is mapped
                    write_reads.append(read)
                    # check if write_reads is too large
                    if len(write_reads) > 1000000:
                        for read in write_reads:
                            output_sam.write(read)
                        write_reads = []

            for read in write_reads:
                output_sam.write(read)


def split_xa_tags(input_sam, output_sam):
    # Open the input SAM file
    write_reads = []
    with pysam.AlignmentFile(input_sam, "r") as infile, \
            pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:

        for read in infile.fetch(until_eof=True):
            # Extract the XA tag if it exists
            xa_tag = read.get_tag('XA') if read.has_tag('XA') else None

            # Write the primary alignment (without the XA tag)
            if xa_tag:
                read.set_tag('XA', None, value_type='Z')  # Remove XA tag
            write_reads.append(read)

            # If the XA tag exists, split and create new alignments
            if xa_tag and not isinstance(xa_tag, int):

                alt_alignments = xa_tag.split(';')
                for alt in alt_alignments:
                    if alt:  # Non-empty alternative alignment
                        # Parse the alternative alignment
                        split = alt.split(',')
                        if len(split) != 4:
                            print(f'Invalid XA tag: {alt}')
                            exit(1)
                        ref_name, pos, cigar, nm = alt.split(',')
                        pos = int(pos)
                        nm = int(nm)

                        # Create a new AlignedSegment object
                        new_read = pysam.AlignedSegment()
                        new_read.query_name = read.query_name  # Same query name

                        new_read.flag = 0
                        # Reference index from reference name
                        new_read.reference_id = infile.get_tid(ref_name)
                        new_read.reference_start = abs(
                            pos) - 1  # Adjust to 0-based start
                        new_read.mapping_quality = read.mapping_quality  # Same MAPQ
                        new_read.cigarstring = cigar  # New CIGAR from XA tag
                        new_read.next_reference_id = read.next_reference_id  # Same RNEXT
                        new_read.next_reference_start = read.next_reference_start  # Same PNEXT
                        new_read.template_length = read.template_length  # Same TLEN

                        new_read.tags = read.tags  # Copy all tags

                        # Set NM tag for mismatches from the XA tag
                        new_read.set_tag('NM', nm)

                        # Set the 256 supplementary flag for the new alternative alignment
                        new_read.flag |= 256  # Supplementary alignment flag

                        # Check if the alignment is on the reverse strand (indicated by negative pos)
                        if pos < 0:
                            new_read.flag |= 16  # Set the reverse complement flag

                        # Write the new alternative alignment to the output
                        write_reads.append(new_read)
            # check if write_reads is too large
            if len(write_reads) > 1000000:
                for read in write_reads:
                    outfile.write(read)
                write_reads = []

        for read in write_reads:
            outfile.write(read)


def getNM(read):
    '''
    Get the NM tag from the read
    '''
    try:
        return read.get_tag("NM")
    except KeyError:
        print("No NM tag found for read: ", read.query_name)
        exit(1)


def keep_only_best_NM(input, output):
    '''
    for each query keep only alignments with the best NM
    '''
    samfile = pysam.AlignmentFile(input, "r")
    current_query = None
    current_best_NM = None
    current_best_alignments = []

    output_file = pysam.AlignmentFile(output, "w", template=samfile)

    for read in samfile:

        if read.query_name != current_query:
            if current_query is not None:
                for alignment in current_best_alignments:
                    output_file.write(alignment)
            current_query = read.query_name
            current_best_NM = getNM(read)
            current_best_alignments = [read]
        else:
            NM = getNM(read)
            if NM < current_best_NM:
                current_best_NM = NM
                current_best_alignments = [read]
            elif NM == current_best_NM:
                current_best_alignments.append(read)

    # write the last query
    for alignment in current_best_alignments:
        output_file.write(alignment)


def get_all_alignments_next_query(samfile, nextLine):
    """
    Retrieve all alignments for the next query name in the given SAM file.

    Args:
        samfile: A pysam.AlignmentFile object representing the SAM file.
        nextLine: The next line in the SAM file to be processed. This line should be the first line of the next query name.

    Returns:
            - reads (list): A list of reads associated with the current query name.
    """
    if nextLine is None:
        current_query_name = None
        current_reads = []
    else:
        current_query_name = nextLine.query_name
        current_reads = [nextLine]

    for read in samfile:
        # Include all reads, regardless of whether they are mapped or unmapped
        if current_query_name is None:
            current_query_name = read.query_name
            current_reads.append(read)
        elif read.query_name == current_query_name:
            # Continue adding reads with the same query name
            current_reads.append(read)
        else:
            # We have reached a new query name; return the current results
            return current_reads, read

    # If we exit the loop and have collected reads, return them
    if current_query_name is not None:
        return current_reads, None

    # Return None if no reads were found
    return [], None


def clean(input_file, output_dir=None):
    # check if cleaned version exists
    if output_dir is not None:
        base = os.path.basename(input_file)
        cleaned = os.path.join(output_dir, os.path.basename(input_file))
        # remove .sam from clened and add _cleaned.sam
        cleaned = os.path.splitext(cleaned)[0] + "_cleaned.sam"
        print(output_dir)
    else:
        base = os.path.splitext(input_file)[0]  # This removes the .sam part
        cleaned = base + "_cleaned.sam"

    print("Cleaned file: ", cleaned)

    if os.path.exists(cleaned):
        print("Cleaned file already exists: ", cleaned)
        return cleaned
    update_cigar_names = ["yara", "rlc"]
    update_cigar_bool = any(name in base.lower()
                            for name in update_cigar_names)

    remove_duplicate_names = ["razer"]
    remove_duplicates_bool = any(name in base.lower()
                                 for name in remove_duplicate_names)

    # Generate a unique identifier based on current time and a random number
    unique_id = f"{int(time.time())}_{random.randint(1000, 9999)}"

    # Create filenames using the unique identifier
    temp1 = f"temp1_{unique_id}.sam"
    temp2 = f"temp2_{unique_id}.sam"

    file_for_sort = temp1 if update_cigar_bool else input_file

    if update_cigar_bool:
        print("Updating CIGAR strings if none given...")
        update_cigar(input_file, temp1)

    print("Sorting file...")
    sort_by_query(file_for_sort, temp2)

    if remove_duplicates_bool:
        print("Removing duplicates...")
        remove_duplicates(temp2, temp1)
    else:
        temp1, temp2 = temp2, temp1

    print("Removing unmapped reads")
    remove_unmapped(temp1, temp2)

    print("Splitting XA tags")
    split_xa_tags(temp2, temp1)

    print("Keeping only best NM alignments")
    keep_only_best_NM(temp1, cleaned)

    os.remove(temp1)
    os.remove(temp2)
    return cleaned


def write_double_value(f, key, value1, value2, file1_name, file2_name):
    f.write(key + "\n")
    f.write(file1_name + ": " + str(len(value1)) + " alignments\n")
    for read in value1:
        f.write(read.to_string() + "\n")
    f.write(file2_name + ": " + str(len(value2)) + " alignments\n")
    for read in value2:
        f.write(read.to_string() + "\n")
    f.write("\n")


def write_doubled_value_dict(dictionary, file1_name, file2_name, output_file_name):
    with open(output_file_name, "w") as output_file:
        for key, value in dictionary.items():
            write_double_value(
                output_file, key, value[0], value[1], file1_name, file2_name)


def write_single_value_dict(dictionary, output_file_name):
    with open(output_file_name, "w") as output_file:
        for _, value in dictionary.items():
            for read in value:
                output_file.write(read.to_string() + "\n")


def find_difference(input1, input2, output_dir=None):
    '''Find the reads that do not have the same best NM alignments in both files'''
    clean1 = clean(input1, output_dir)
    clean2 = clean(input2, output_dir)
    print("Cleaned files: ", clean1, clean2)

    best1 = pysam.AlignmentFile(clean1, "r")
    best2 = pysam.AlignmentFile(clean2, "r")

    base1 = os.path.basename(input1)
    base2 = os.path.basename(input2)

    # find the directory in whcih input1 and input2 are located
    dir_name = os.path.dirname(input1)
    if output_dir is None:
        output_dir = os.path.join(
             dir_name, "compare_" + base1 + "_" + base2)

    # make directory if not exists
    if not os.path.exists(output_dir):
        print("Creating output directory: ", output_dir)
        os.makedirs(output_dir)
    else:
        print("Output directory already exists: ", output_dir)
        # read the meta file for the two numbers
        meta_file = os.path.join(output_dir, "meta.txt")
        number_file1_more, number_file2_more = None, None
        if os.path.exists(meta_file):
            with open(meta_file, "r") as f:
                lines = f.readlines()
                # check lines 2 and 4 for the two numbers
                number_file1_more = int(lines[1].strip())
                number_file2_more = int(lines[3].strip())

        return output_dir, number_file1_more, number_file2_more

    print("Output directory: ", output_dir)

    # Prepare output files dynamically
    output_only_file1 = open(os.path.join(output_dir, "only_" + base1), "w")
    output_only_file2 = open(os.path.join(output_dir, "only_" + base2), "w")
    output_score1_better = open(os.path.join(
        output_dir, "score_better_" + base1), "w")
    output_score2_better = open(os.path.join(
        output_dir, "score_better_" + base2), "w")
    output_file1_more = open(os.path.join(
        output_dir, "more_alignments_" + base1), "w")
    output_file2_more = open(os.path.join(
        output_dir, "more_alignments_" + base2), "w")

    output_file_not_same = open(os.path.join(
        output_dir, "not_same_" + base1 + "_" + base2), "w")

    iteration = 0

    nextLine1 = None
    nextLine2 = None

    mapped1 = 0
    mapped2 = 0

    advance1 = True
    advance2 = True

    number_of_reads = 0

    number_file1_more = 0
    number_file2_more = 0

    while True:

        print(
            f"\rReads processed in file 1: {number_of_reads}", end='', flush=True)
        iteration += 1

        if advance1:
            reads1, nextLine1 = get_all_alignments_next_query(best1, nextLine1)
            mapped1 += 1 if len(reads1) > 0 else 0
            number_of_reads += 1
        if advance2:
            reads2, nextLine2 = get_all_alignments_next_query(best2, nextLine2)
            mapped2 += 1 if len(reads2) > 0 else 0

        if len(reads1) == 0 and len(reads2) == 0:
            break

        # default: advance both
        advance1 = True
        advance2 = True

        # Handle read comparison
        if len(reads1) == 0:
            # no more alignments in file1, so all alignments unique to file2
            for read in reads2:
                output_only_file2.write(read.to_string() + "\n")
            advance1 = False  # do not advance file1
            continue
        if len(reads2) == 0:
            # no more alignments in file2, so all alignments unique to file1
            for read in reads1:
                output_only_file1.write(read.to_string() + "\n")
            advance2 = False  # do not advance file2
            continue

        if reads1[0].query_name == reads2[0].query_name:
            # same query name, read aligned in both files

            # get the edit distance from the NM tag
            score1 = reads1[0].get_tag("NM")
            # get the edit distance from the NM tag
            score2 = reads2[0].get_tag("NM")
            if score1 < score2:
                # different NM scores... file1 is better
                write_double_value(
                    output_score1_better, reads1[0].query_name, reads1, reads2, base1, base2)
            elif score1 > score2:
                # different NM scores... file2 is better
                write_double_value(
                    output_score2_better, reads1[0].query_name, reads1, reads2, base1, base2)
            elif len(reads1) < len(reads2):
                # same NM score, but file1 has fewer alignments
                number_file2_more += 1
                write_double_value(
                    output_file2_more, reads1[0].query_name, reads1, reads2, base1, base2)
            elif len(reads1) > len(reads2):
                # same NM score, but file1 has more alignments
                number_file1_more += 1
                write_double_value(
                    output_file1_more, reads1[0].query_name, reads1, reads2, base1, base2)
            else:
                # same number of alignments and same NM score for all alignments
                # we need to check if the alignments are rougly at the same position
                ref_pos_dict1 = {}
                ref_pos_dict2 = {}
                for read in reads1:
                    ref_id = read.reference_name
                    if ref_id not in ref_pos_dict1:
                        ref_pos_dict1[ref_id] = []
                    ref_pos_dict1[ref_id].append(read.reference_start)
                for read in reads2:
                    ref_id = read.reference_name
                    if ref_id not in ref_pos_dict2:
                        ref_pos_dict2[ref_id] = []
                    ref_pos_dict2[ref_id].append(read.reference_start)

                # check if the reference positions are roughly the same
                if ref_pos_dict1.keys() != ref_pos_dict2.keys():
                    # different reference IDs, so we cannot compare
                    write_double_value(
                        output_file_not_same, reads1[0].query_name, reads1, reads2, base1, base2)
                else:
                    for ref_id in ref_pos_dict1.keys():
                        pos1 = ref_pos_dict1[ref_id]
                        pos2 = ref_pos_dict2[ref_id]
                        if len(pos1) != len(pos2):
                            # different number of alignments for this reference ID, so we cannot compare
                            write_double_value(
                                output_file_not_same, reads1[0].query_name, reads1, reads2, base1, base2)
                            break
                        # check if the positions are roughly the same
                        pos1.sort()
                        pos2.sort()
                        for p1, p2 in zip(pos1, pos2):
                            if abs(p1 - p2) > 50:
                                # positions are not roughly the same so write to output
                                write_double_value(
                                    output_file_not_same, reads1[0].query_name, reads1, reads2, base1, base2)
                                break
            continue

        # different query names, one of the files has alignments for a query while the other does not
        # due to the sorting by query name, we can compare the first read of each file
        # and determine which query is missing in which file

        # strip the last part of the query name and convert to integer
        query1 = int(reads1[0].query_name.split(".")[-1])
        query2 = int(reads2[0].query_name.split(".")[-1])

        if query1 < query2:
            # file 1 has an alignment for query1, this query is not present in file2
            for read in reads1:
                output_only_file1.write(read.to_string() + "\n")
            advance2 = False  # do not advance file2
            continue

        if query1 > query2:
            # file 2 has an alignment for query2, this query is not present in file1
            for read in reads2:
                output_only_file2.write(read.to_string() + "\n")
            advance1 = False  # do not advance file1
            continue
    print()  # Print a newline after the progress output
    print(f"Mapped reads in {base1}: ", mapped1)
    print(f"Mapped reads in {base2}: ", mapped2)

    # Close output files
    output_only_file1.close()
    output_only_file2.close()
    output_score1_better.close()
    output_score2_better.close()
    output_file1_more.close()
    output_file2_more.close()

    # write the two numbers to a meta file
    meta_file = os.path.join(output_dir, "meta.txt")
    with open(meta_file, "w") as f:
        f.write("Number more in file 1:\n" + str(number_file1_more) + "\n")
        f.write("Number more in file 2:\n" + str(number_file2_more) + "\n")

    return output_dir, number_file1_more, number_file2_more


def printUsage():
    print(
        "Usage: compare_sam.py <file1.sam> <file2.sam> <reference_fasta> [<max_edit_distance> <max_error_rate>]")


def printline():
    print()
    print("-" * 80)
    print()


def print_double_line():
    print()
    print("=" * 80)
    print()


def get_tool(base):
    '''
    Get the tool name from the base name of the file
    '''
    if "yara" in base.lower():
        return "YARA"
    elif "rlc" in base.lower():
        return "Columba RLC"
    elif "razer" in base.lower():
        return "Razer"
    elif "columba" in base.lower() or "vanilla" in base.lower():
        return "Columba"
    elif "mem" in base.lower():
        return "BWA-mem"
    else:
        return base


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4 or len(sys.argv) > 7:
        printUsage()
        sys.exit(1)
    file1_name = sys.argv[1]
    file2_name = sys.argv[2]

    fasta = sys.argv[3]

    edit_distance = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    if edit_distance < 0:
        print("Edit distance must be a non-negative integer.")
        sys.exit(1)
    threshold = float(sys.argv[5]) if len(sys.argv) > 5 else 0.08

    number_reads = int(sys.argv[6]) if len(sys.argv) > 6 else pow(10, 6)
    number_reads_same = number_reads

    tolerance = 2 * edit_distance + 1

    printline()
    print(f"Comparing {file1_name} and {file2_name}...")

    output_dir, number_file1_more, number_file2_more = find_difference(
        sys.argv[1], sys.argv[2])

    print(f"{output_dir} contains the results of the comparison between {file1_name} and {file2_name}.")

    print("Inpsection of these files will now be performed. You can always inspect the files manually."
          " If you want to skip this step, just press Ctrl+C.")

    printline()

    tools = [get_tool(file1_name), get_tool(file2_name)]

    print("Loading sequences from the FASTA file for non-ACGT character inspection...")

    sequences = load_sequences(fasta)

    # load header_obj from file1_name
    header_obj = pysam.AlignmentFile(file1_name, "r").header

    files = os.listdir(output_dir)

    # remove files that end with .txt, as these are not relevant for the inspection
    files = [f for f in files if not (
        f.endswith(".txt") or f.endswith(".png"))]

    print_double_line()
    # step 1: inspect the "more_alignments" files
    print("Inspecting 'more_alignments' files for reasons for why one tool has more alignments for a given read than the other...\n")
    files_more = [f for f in files if ((f.startswith("more_alignments_")) and (
        ("_distributions" or ".png") not in f))]

    for filename in files_more:
        number_more = number_file1_more
        tool1, tool2 = tools[0], tools[1]
        # check if filename contains base file 1
        if not os.path.basename(file1_name) in filename:
            number_more = number_file2_more
            tool1, tool2 = tool2, tool1
        print(
            f"Inspecting alignments for reads with more alignments reported by {tool1} than in {tool2}  ({filename})...")
        numbers_dict, amounts = inspect_more_alignments.parse_file(os.path.join(
            output_dir, filename), header_obj, tolerance, sequences, number_more)
        print(f"\rResults for {filename}:")
        print(f"  Total reads inspected: {numbers_dict['count']}")
        print(
            f"  Alignments by {tool2} are subset of alignments by {tool1}: {numbers_dict['subset']}")
        print(
            f"    ├─ Explained by non-ACGT characters: {numbers_dict['subset_all_non_acgt']}")
        print(
            f"    ├─ Explained by redundancy handling: {numbers_dict['subset_all_close']}")
        print(
            f"    ├─ Explained by error rate calculation: {numbers_dict['subset_error_rate']}")
        print(
            f"    ├─ Explained by short reads: {numbers_dict['subset_short_read']}")
        print(
            f"    ├─ Explained by CIGAR clipping: {numbers_dict['subset_all_s_or_h']}")
        print(f"    └─ Unexplained alignments: {numbers_dict['unexplained']}")
        if (numbers_dict["count"] != numbers_dict["subset"]):
            print(
                f"  Alignments by {tool2} are not a subset of alignments by {tool1}: {numbers_dict['count'] - numbers_dict['subset']}")
            print(
                f"    └─  Explained by non-ACGT characters: {numbers_dict['non_subset_all_non_acgt']}")

        # subtract the number of reads inspected here
        number_reads_same -= numbers_dict['count']

        inspect_more_alignments.plot_distributions(
            amounts, tool1, tool2, os.path.join(output_dir, filename + "_distributions.png"))
        # chech if last in the loop
        if filename == files_more[-1]:
            print_double_line()
        else:
            printline()

    # step 2: inspect the "not_same" files
    # As these have always been empty for the tests, we do not have automation for this step.
    # some of the steps will probably be the same as in more_alignments and/or score better
    # not_same implies the same number of alignments with same best NM, but different reference positions
    # which means both tools report an alignment that the other tool does not report

    # report number of times the word alignments appears in the file and report the half of this number
    not_same_file = os.path.join(
        output_dir, "not_same_" + os.path.basename(file1_name) + "_" + os.path.basename(file2_name))
    if os.path.exists(not_same_file):
        print("Inspecting 'not_same' file, for reads that have the same best NM tag and the same number of alignments, but different reference positions (highly unlikely scenario)...")
        with open(not_same_file, "r") as f:
            lines = f.readlines()
            alignments_count = sum(1 for line in lines if "alignments" in line)
            print(
                f"Number of reads in {not_same_file}: {alignments_count // 2}")
            # subtract the number of reads inspected here
            number_reads_same -= alignments_count // 2

    print_double_line()

    # step 3: inspect the only files
    # these files contain reads that are only aligned in only one of the files
    files_only = [f for f in files if f.startswith("only_")]
    print(f"Inspecting only_ files for reasons why alignments for a given read are only present for one tool...\n")
    for filename in files_only:
        tool = get_tool(os.path.basename(filename))
        other_tool = tools[0] if tool == tools[1] else tools[1]
        print(
            f"Inspecting alignments for reads only reported as aligned by {tool} and not by {other_tool}...")

        numbers_dict, unknown_reasons = inspect_only(os.path.join(
            output_dir, filename), header_obj, sequences, threshold)
        print(
            f"Results for alignments for reads only reported as aligned by {tool} and not by {other_tool}:")
        print(f"  Total reads inspected: {numbers_dict['count']}")
        print(
            f"    ├─ Reads with S/H in CIGAR (clipping): {numbers_dict['SorH']}")
        print(
            f"    ├─ Reads with non-ACGT characters: {numbers_dict['non_acgt_characters']}")
        print(
            f"    ├─ Reads exceeding error rate: {numbers_dict['error_rate_above_threshold']}")
        print(
            f"    ├─ Reads with short read length: {numbers_dict['short read']}")
        print(f"    └─ Reads with unknown reasons: {len(unknown_reasons)}")
        unknown_file = open(os.path.join(
            output_dir, f"{filename}_unknown_reasons.txt"), "w")
        # write unknown reasons to file
        for read_id, alignments in unknown_reasons:
            for aln in alignments:
                unknown_file.write(f"{aln.query_name}\n")
        unknown_file.close()
        print(
            f"Query ids for unknown reason missing alignments are written to {unknown_file.name}")

        # subtract the number of reads inspected here
        number_reads_same -= numbers_dict['count']

        if filename == files_only[-1]:
            print_double_line()
        else:
            printline()

    # step 4: inspect the score better files
    # these files contain reads that have a better NM score in one of the files
    files_score_better = [f for f in files if f.startswith("score_better_")]
    print(f"Inspecting score_better_ files for edit distance score differences...\n")
    for filename in files_score_better:

        tool = get_tool(os.path.basename(filename))
        other_tool = tools[0] if tool == tools[1] else tools[1]
        print(tool, other_tool, tools)
        print(
            f"Inspecting reads with alignments with better edit distance score reported by {tool} than by {other_tool}...")
        unknown_reasons, relative_same_positions, numbers_dict = inspect_better_score.parse_file(
            os.path.join(output_dir, filename), header_obj, tolerance, sequences)

        print(
            f"Results for reads with alignments with better edit distance socre reported by {tool} than by {other_tool}:")
        print(f"  Total reads inspected: {numbers_dict['count']}")
        print(
            f"    ├─ Reads with S/H in CIGAR (clipping): {numbers_dict['SorH']}")
        print(
            f"    ├─ Relative same positions (within tolerance): {numbers_dict['relative_same_positions']}")
        print(
            f"    │   └─ With non-ACGT characters: {numbers_dict['non_acgt_handling_rel_same_pos']}")
        print(
            f"    └─ Widely different alignments (missed optimal): {numbers_dict['widely_different_alignment']}")
        print(
            f"        └─ With non-ACGT characters: {numbers_dict['non_acgt_handling_widely_different']}\n")

        # subtract the number of reads inspected here
        number_reads_same -= numbers_dict['count']

        file_score_better_ids = open(os.path.join(
            output_dir, "score_better_ids.txt"), "w")
        print("Writing IDs of reads with better score not explained by non-ACGT or clipping in CIGAR string to file...")
        file_score_better_ids.write(
            "Relative same positions but no non-ACGT character in reference:\n")
        for read_id, a1, a2 in relative_same_positions:
            file_score_better_ids.write(f"Read ID: {read_id}")
            file_score_better_ids.write("Tool1:")
            for seq_id, pos, cigar in a1:
                file_score_better_ids.write(f"\t{seq_id} {pos} {cigar}")
            file_score_better_ids.write("Tool2:")
            for seq_id, pos, cigar in a2:
                file_score_better_ids.write(f"\t{seq_id} {pos} {cigar}")

        file_score_better_ids.write("Unknown reasons (no non-ACGT):\n")
        for read_id, a1, a2 in unknown_reasons:
            file_score_better_ids.write(f"Read ID: {read_id}")
            file_score_better_ids.write("Tool1:")
            for seq_id, pos, cigar in a1:
                file_score_better_ids.write(f"\t{seq_id} {pos} {cigar}")
            file_score_better_ids.write("Tool2:")
            for seq_id, pos, cigar in a2:
                file_score_better_ids.write(f"\t{seq_id} {pos} {cigar}")

        if filename == files_score_better[-1]:
            print_double_line()
        else:
            printline()

    print(
        f"Number of reads with same alignments (+- {tolerance} bp) in both files: {number_reads_same}")
