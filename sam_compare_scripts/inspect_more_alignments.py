"""
This script analyzes files named more_alignments_tool1_tool2.sam, which contain alignment data
for reads processed by two different alignment tools.

Input file structure per read is:

    Read ID line

    Info line for tool 1 indicating number of alignments

    Alignment lines for tool 1 (X lines)

    Info line for tool 2 indicating number of alignments

    Alignment lines for tool 2 (Y lines)

    Blank line separating entries

The script parses these files, compares the sets of alignments reported by each tool for each read,
and categorizes differences based on sequence content (e.g., presence of non-ACGT characters),
positional similarity, alignment error rates, read length, and presence of soft/hard clipping in CIGAR strings.

It processes the input in batches using parallel workers to improve performance, using pysam to
parse SAM alignments and leveraging global shared data to avoid redundant serialization.

The output includes statistics on how often the alignments from tool 2 are subsets of tool 1's,
along with explanations for differences where possible.

It is intended as a helper for compare_sam.py and assumes that compare_sam.py's initial processing
step has already been performed. Generally, it does not need to be run separately.

Usage:
python find_missing_alignments.py <input_directory> <sam_file_with_header> <edit_distance>
"""

import pysam
import os
import sys
import subprocess
import matplotlib.pyplot as plt


from concurrent.futures import as_completed, ProcessPoolExecutor


from utils import (
    contains_non_acgt,
    make_dict,
    match_other_as_dict,
    alignment_length,
    read_length,
    s_or_h_in_cigar
)

def match(pos, lst, tolerance):
    return any(r == pos[0] and abs(rp - pos[1]) <= tolerance for r, rp, _, _ in lst)


def match_parallel(alignments_2, alignments1_dict, tolerance):
    return [match_other_as_dict(a, alignments1_dict, tolerance) for a in alignments_2]


def different_error_rate(cigar_tuples, edit_distance):
    """
    Check if calculation of error rate via alignment length differs from calculation via read length
    """
    a_length = alignment_length(cigar_tuples)
    r_length = read_length(cigar_tuples)

   

    error_rate_alignment = (edit_distance / a_length) * 100
    error_rate_read = (edit_distance / r_length) * 100
    # floor both values to integer and see if they differ
    return int(error_rate_alignment) != int(error_rate_read)


def check_reason(
    read_id, distinct_alignments, dict2, tolerance, sequences, lossless, amount_more
):
    """
    Process a chunk of read IDs to find distinct alignments and categorize them.
    """

    aln_format = [(a[0], a[1], a[2]) for a in distinct_alignments]
    cat = ""

    if contains_non_acgt(aln_format, sequences):
        cat = "subset_all_non_acgt"
    elif all(match_other_as_dict(a, dict2, tolerance) for a in distinct_alignments):
        cat = "subset_all_close"
    elif lossless and all(
        different_error_rate(a[2], a[3]) for a in distinct_alignments
    ):
        cat = "subset_error_rate"
    elif all(read_length(a[2]) < 100 for a in distinct_alignments):
        cat = "subset_short_read"
    elif all(s_or_h_in_cigar(a[2]) for a in distinct_alignments):
        cat = "subset_all_s_or_h"
    else:
        cat = "unexplained"

    return cat, amount_more


def read_blocks(f):
    """
    Generator to read blocks of alignments from the file.
    Yields tuples of (read_id, alignments_1, alignments_2).
    The alignments are given as lists of strings.
    Each block consists of:
    1. Read ID
    2. Info line for tool 1 with number of alignments
    3. Alignment lines for tool 1 (probably multiple lines)
    4. Info line for tool 2 with number of alignments
    5. Alignment lines for tool 2 (probably multiple lines)
    6. Blank line
    """
    current_id = None
    current_n_1 = 0
    current_n_2 = 0
    current_aln_1 = []
    current_aln_2 = []

    i =0
    for line in f:
        line = line.rstrip('\n')

        if current_id is None:
            current_id = line

        elif current_n_1 == 0:
            current_n_1 = int(line.split(":")[1].split()[0])

        elif len(current_aln_1) < current_n_1:
            current_aln_1.append(line)

        elif current_n_2 == 0:
            current_n_2 = int(line.split(":")[1].split()[0])

        elif len(current_aln_2) < current_n_2:
            current_aln_2.append(line)

        else:
            # assert that line is blank
            if line.strip() != "":
                raise ValueError(f"Expected blank line, got: {line}")
            # Expect blank line here - yield block
            yield (current_id, current_aln_1, current_aln_2)

            

            # Reset state for next block
            current_id = None
            current_n_1 = 0
            current_n_2 = 0
            current_aln_1 = []
            current_aln_2 = []

            

            
        

    # Yield last block if file ends mid-block
    if current_id is not None:
        yield (current_id, current_aln_1, current_aln_2)


# --- Single read processing (now works on plain tuples) ---
def process_one_read(args, tolerance, sequences, lossless):
    read_id, alignments_1, alignments_2 = args

    if len(alignments_1) < len(alignments_2):
        alignments_1, alignments_2 = alignments_2, alignments_1

    dict2 = make_dict(alignments_2)
    match_bool = match_parallel(alignments_1, dict2, tolerance)
    distinct_alignments = [a for a, m in zip(
        alignments_1, match_bool) if not m]
    amount_more = len(alignments_1) - len(alignments_2)

    dict1 = make_dict(alignments_1)
    subset = all(match_parallel(alignments_2, dict1, tolerance))

    if subset:
        cat, amount_more = check_reason(
            read_id,
            distinct_alignments,
            dict2,
            tolerance,
            sequences,
            lossless,
            amount_more,
        )
    else:
        cat = "unmatched"
        aln_format = [
            (key, pos, cigartuples)
            for key, values in dict2.items()
            for pos, cigartuples, _ in values
        ]
        if contains_non_acgt(aln_format, sequences):
            cat = "non_subset_all_non_acgt"
        
        
    return subset, cat, amount_more


def lines_to_tuples(aln_raw, header):
    """
    Convert raw alignment lines to tuples of (reference_name, reference_start, cigartuples, edit_distance).
    Each line is parsed using the provided header to create a pysam.AlignedSegment object.
    This function assumes that aln_raw is an iterable of alignment lines.
    Each line should be a valid SAM/BAM format string.
    The header is a pysam.AlignmentHeader object.
    """
    return [
        (
            a.reference_name,
            a.reference_start,
            a.cigartuples,
            a.get_tag("NM")
        )
        for a in (pysam.AlignedSegment.fromstring(line.strip(), header) for line in aln_raw)
    ]


def convert_strings_to_tuples(batch_raw,  header):
    """
    Convert a raw batch containing (read_id, aln_1_string, aln_2_string) tuples to a batch containing (read_id, aln_1_tuples, aln_2_tuples).
    Uses lines_to_tuples to convert the alignment strings to tuples.
    """
    return [
        (read_id, lines_to_tuples(aln_1, header), lines_to_tuples(aln_2, header))
        for read_id, aln_1, aln_2 in batch_raw
    ]


_global_sequences = None
_global_header = None

import traceback
def init_worker(sequences, header_dict):
    """
    Initialize global variables for each worker process.
    This avoids repeatedly passing large data structures.
    """
    global _global_sequences, _global_header
    _global_sequences = sequences
    _global_header = pysam.AlignmentHeader.from_dict(header_dict)


def process_batch(args_batch_raw, tolerance, lossless, index, n_batches):
    results = []
    global _global_sequences, _global_header
  
    try:
        
        args_batch = convert_strings_to_tuples(args_batch_raw, _global_header)
        
        for args in args_batch:
            result = process_one_read(
                args, tolerance, _global_sequences, lossless)
            subset, cat, amount_more = result
            results.append((subset, cat, amount_more))

        print(f"\rProcessed batch {index}/{n_batches} with {len(args_batch)} reads", end="", flush=True)

        assert len(results) == len(args_batch), "Batch processing mismatch"
        return results, index

    except Exception as e:
        print("exception in batch processing:", e)
        with open(f"error_{index}.log", "w") as f:
            traceback.print_exc(file=f)
        raise RuntimeError(f"Error in batch {index}") 


def submit_batch(executor, batch_raw, tolerance, lossless, batch_index, futures, n_batches):

    if len(batch_raw) > 0:
        future = executor.submit(
            process_batch, batch_raw, tolerance, lossless, batch_index, n_batches)
        futures.add(future)


def process_futures(futures, numbers_dict, amounts_per_category, total_reads):
    """
    Collect and process results from completed futures.
    Updates counts and accumulates statistics.
    """
    for future in as_completed(futures):
        
        try:
            result_batch, index = future.result()
          
           
            for result_read in result_batch:
                subset, cat, amount_more = result_read
                if subset:
                    numbers_dict["subset"] += 1
                numbers_dict[cat] += 1
                amounts_per_category[cat].append(amount_more)
                numbers_dict["count"] += 1

                if numbers_dict["count"] % 16 == 0:
                    print(
                        f"\rProcessed ({numbers_dict['count']}/{total_reads})", end="", flush=True)
        except Exception as e:
            print(f"Worker error: {e}")


def parse_file(path, header, tolerance, sequences, total_reads=None):
    """
    Main function to parse the alignment file and process reads in parallel.

    Workflow:
    1. Determine total reads from file using grep.
    2. Open the file and read in blocks using a generator.
    3. Batch blocks and submit to ProcessPoolExecutor workers.
    4. Use a global sequences object to avoid repeated serialization.
    5. Collect results and update summary statistics.
    """
    numbers_dict = {
        "subset": 0,
        "subset_all_non_acgt": 0,
        "subset_all_close": 0,
        "subset_error_rate": 0,
        "subset_short_read": 0,
        "unknown_reason": 0,
        "unexplained": 0,
        "unmatched": 0,  # This is for the case where alignments of tool 2 are not a subset of tool 1
        "count": 0,
        "subset_all_s_or_h": 0,  
        "non_subset_all_non_acgt": 0,  # This is for the case where alignments of tool 2 are not a subset of tool 1, but contain non-ACGT characters
    }
    amounts_per_category = {key: [] for key in numbers_dict}

    lossless = "bwa_mem" not in path and "bowtie2" not in path

    # Count total reads (assuming 2 alignments per read block)
    if total_reads is None:
        # Use grep to count the number of alignments in the file
        # This assumes that each alignment line contains the word "alignments"
        # and that there are two alignments per read (one for each tool).
        print(f"Counting total reads in {path}...")
        result = subprocess.run(
            ["grep", "-c", "alignments", path], capture_output=True, text=True)
        total_reads = int(result.stdout.strip()) // 2
    print(f"Total reads to process: {total_reads}")

    if total_reads == 0:
        print("No reads found.")
        return numbers_dict, amounts_per_category

    BATCH_SIZE = 128
    max_workers = 32
    n_batches = total_reads // BATCH_SIZE + 1

    futures = set()

    header_dict = header.to_dict()

    # Create a manager and a shared counter
 


    with open(path) as f, ProcessPoolExecutor(
        max_workers=max_workers, initializer=init_worker, initargs=(
            sequences, header_dict,)
    ) as executor:
        args_batch_raw = []
        batch_index = 0

        for block in read_blocks(f):
            args_batch_raw.append(block)

            if len(args_batch_raw) == BATCH_SIZE:
            

                submit_batch(executor, args_batch_raw[:], tolerance, lossless,
                             batch_index, futures, n_batches)
                batch_index += 1
                args_batch_raw.clear()

        # Submit remaining reads if any
        if args_batch_raw:
             submit_batch(executor, args_batch_raw, tolerance, lossless,
                             batch_index, futures, n_batches)

        # Wait for all remaining futures
        process_futures(futures, numbers_dict,
                        amounts_per_category, total_reads)

    return numbers_dict, amounts_per_category


def plot_distributions(amounts_per_category, tool1, tool2, destination=None):
    """
    Plot the distributions of the amounts per category.
    """



    for category, amounts in amounts_per_category.items():

        if not amounts:
            continue  # Skip categories with no amounts
        plt.hist(amounts, bins=20, alpha=0.5, label=category)

    plt.xlabel(f"Number of missed alignments per category by {tool2}")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of missed alignments per category by {tool2} that were reported by {tool1}")
    plt.legend()
    if destination:
        plt.savefig(destination)
        return
    plt.tight_layout()
    plt.show()


# Usage
if __name__ == "__main__":

    # get input directory from user
    if len(sys.argv) != 4:
        print(
            "Usage: python find_missing_alignments.py <input_directory> <sam_file_with_header> <edit distance>"
        )
        sys.exit(1)

    input_directory = sys.argv[1]
    sam_file_with_header = sys.argv[2]

    ref_sam = pysam.AlignmentFile(sam_file_with_header, "r")  # or "rb" for BAM
    header_obj = ref_sam.header  # this is already a pysam.AlignmentHeader
    edit_distance = int(sys.argv[3])

    tolerance = 2 * edit_distance + 1  # tolerance for matching positions

    # do this for all files in directory starting with "more_alignments_"

    for filename in os.listdir(input_directory):

        if filename.startswith("more_alignments_"):
            numbers_dict, amounts_per_category = parse_file(
                os.path.join(input_directory, filename), header_obj, tolerance
            )
            plot_distributions(amounts_per_category)

            print(f"\rResults for {filename}:")
            print(f"  Total reads inspected: {numbers_dict['count']}")
            print(
                f"  Alignments of tool 2 are subset of alignments in tool 1: {numbers_dict['subset']}"
            )
            print(
                f"    ├─ Explained by non-ACGT characters: {numbers_dict['subset_all_non_acgt']}"
            )
            print(
                f"    ├─ Explained by redundancy handling: {numbers_dict['subset_all_close']}"
            )
            print(
                f"    ├─ Explained by error rate calculation: {numbers_dict['subset_error_rate']}"
            )
            print(
                f"    ├─ Explained by short reads: {numbers_dict['subset_short_read']}"
            )
            print(
                f"    └─ Unexplained alignments: {numbers_dict['unmatched']}")

            input("Press Enter to continue...")
