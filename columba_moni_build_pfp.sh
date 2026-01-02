#!/bin/bash

# Script: columba_build_pfp.sh
# Description: This script performs the Columba build process for PFP.
# Author: Lore Depuydt - lore.depuydt@ugent.be

# -----------------------------------------------------------------------------
# INITIALIZATION & DEFAULTS
# -----------------------------------------------------------------------------

# Executable paths
columba_build_exe="./columba_build"
big_bwt_exe="./../external/Big-BWT/bigbwt"
moni_align_exe="./../external/moni-align/build"

# Default seed length
seedLength=100

# Optional Big-BWT parameters
ws=0 
mod=0

# Moni-Align specific variables
moni_ref=""
moni_vcf=""
moni_samples=""
moni_output=""

# Array to store fasta files
fasta_files=()

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

showUsage() {
    echo "Usage: $0 [options] -r <index_name>"
    echo
    echo "Required arguments (Columba):"
    echo "  -r <index_name>      Name/location of the COLUMBA index to be created."
    echo
    echo "Optional arguments (Columba):"
    echo "  -f <fasta_files>     Space-separated list of FASTA files."
    echo "  -F <fasta_file_list> Path to a file containing a list of FASTA files."
    echo "  -l <seedLength>      Seed length (default: $seedLength)."
    echo "  -w <ws>              Big-BWT window size."
    echo "  -p <mod>             Big-BWT mod value."
    echo
    echo "Moni-Align arguments (Required if running moni build):"
    echo "  -m <ref_fasta>       Reference FASTA for moni-align (-r in moni)."
    echo "  -v <vcf_file>        VCF file for moni-align (-v in moni)."
    echo "  -s <samples_file>    Samples file for moni-align (-S in moni)."
    echo "  -o <moni_out>        Output prefix for moni-align (-o in moni)."
}

runCommandWithTime() {
    local command="$1"
    shift
    (/usr/bin/time -v "$command" "$@") || {
        local status=$?
        echo "Error: Command '$command $@' failed with exit status $status." >&2
        exit $status
    }
}

parseOptions() {
    # Added m, v, s, o to the getopts string
    while getopts ":l:r:f:F:w:p:m:v:s:o:" opt; do
        case $opt in
            l) seedLength=$OPTARG ;;
            r) index_name=$OPTARG ;;
            w) ws=$OPTARG ;;
            p) mod=$OPTARG ;;
            m) moni_ref=$OPTARG ;;
            v) moni_vcf=$OPTARG ;;
            s) moni_samples=$OPTARG ;;
            o) moni_output=$OPTARG ;;
            f)
                fasta_files+=("$OPTARG")
                while [[ $OPTIND -le $# && ! ${!OPTIND} =~ ^- ]]; do
                    fasta_files+=("${!OPTIND}")
                    OPTIND=$((OPTIND + 1))
                done
                ;;
            F)
                if [[ -f $OPTARG ]]; then
                    while IFS= read -r line; do fasta_files+=("$line"); done <"$OPTARG"
                else
                    echo "Error: File '$OPTARG' not found." >&2
                    exit 1
                fi
                ;;
            \?) echo "Invalid option: -$OPTARG" >&2; showUsage; exit 1 ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; showUsage; exit 1 ;;
        esac
    done
    shift $((OPTIND - 1))

    if [ -z "$index_name" ] || [ "${#fasta_files[@]}" -eq 0 ]; then
        showUsage
        exit 1
    fi
}

# -----------------------------------------------------------------------------
# PARSE OPTIONS (Moved Up)
# -----------------------------------------------------------------------------
# We parse options first to determine the index_name for the log folder path.
parseOptions "$@"

# -----------------------------------------------------------------------------
# LOGGING SETUP (Updated)
# -----------------------------------------------------------------------------
# Capture start time
start_time=$(date +%s)
timestamp=$(date +"%Y%m%d_%H%M%S")

# Determine directories based on -r index_name
index_dir=$(dirname "$index_name")
# If index_name has no path, default to current directory
if [ -z "$index_dir" ] || [ "$index_dir" == "." ]; then
    index_dir="."
fi

# Create a specific folder for this run inside the index folder
log_run_dir="${index_dir}/logs/run_${timestamp}"
mkdir -p "$log_run_dir"

# Define the main log file inside that folder
log_filename="${log_run_dir}/columba_moni_build_pfp.log"

echo "-------------------------------------------------------------"
echo "Log file initiated: $log_filename"
echo "All output will be saved to this file."
echo "-------------------------------------------------------------"

# Redirect stdout (1) and stderr (2) to the log file while still showing them on screen
exec > >(tee -a "$log_filename") 2>&1

# -----------------------------------------------------------------------------
# MAIN SCRIPT LOGIC
# -----------------------------------------------------------------------------

echo "Welcome to the Columba build process with prefix-free parsing!"
echo "-------------------------------------------------------------"
echo "Index name: $index_name"
echo "Seed length: $seedLength"
echo "Logs for this run will be stored in: $log_run_dir"

# Execute Moni-Align Build if parameters are provided
# Optimized Moni-Align call
if [ -n "$moni_ref" ] && [ -n "$moni_vcf" ]; then
    echo "Running Moni-Align build..."
    
    # Build the command array dynamically
    moni_cmd=("${moni_align_exe}/moni" build -r "$moni_ref" -v "$moni_vcf" -H12 -o "$moni_output")
    
    # Only add -S if the file is provided
    if [ -n "$moni_samples" ]; then
        moni_cmd+=(-S "$moni_samples")
    fi

    runCommandWithTime "${moni_cmd[@]}"
else
    echo "Skipping Moni-Align (missing -m or -v flags)."
fi
echo "Moni-Align index built!"
echo "-------------------------------------------------------------"


# Preprocessing
echo "Start preprocessing the fasta file(s) with Columba..."
runCommandWithTime "$columba_build_exe" --preprocess -l "$seedLength" -r "$index_name" -f "${fasta_files[@]}"
echo "Preprocessing done!"
echo "-------------------------------------------------------------"

# Big-BWT Setup
base="${index_name}"
big_bwt_args=("$big_bwt_exe" -e -s -v "$base")
[ "$ws" -gt 0 ] && big_bwt_args+=(-w "$ws")
[ "$mod" -gt 0 ] && big_bwt_args+=(-p "$mod")

# Prefix-Free Parsing (Original)
echo "Start prefix-free parsing for the original string..."
runCommandWithTime "${big_bwt_args[@]}"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"


# Prefix-Free Parsing (Reverse)
big_bwt_args=("$big_bwt_exe" -e -s -v "${base}.rev")
[ "$ws" -gt 0 ] && big_bwt_args+=(-w "$ws")
[ "$mod" -gt 0 ] && big_bwt_args+=(-p "$mod")

echo "Start prefix-free parsing for the reverse string..."
runCommandWithTime "${big_bwt_args[@]}"
echo "Prefix-free parsing done (reverse)!"
echo "-------------------------------------------------------------"

# Final Columba Index
echo "Start building the Columba index..."
runCommandWithTime "$columba_build_exe" --pfp -r "$index_name"
echo "Columba index built!"
echo "-------------------------------------------------------------"

# -----------------------------------------------------------------------------
# LOG MANAGEMENT (Subprocess Cleanup)
# -----------------------------------------------------------------------------
echo "Moving subprocess logs to run directory..."

# Function to safely move files if they exist
safe_move_log() {
    local src="$1"
    local dest_dir="$2"
    if [ -f "$src" ]; then
        mv "$src" "$dest_dir"
        echo "Moved log: $src -> $dest_dir"
    else
        echo "Log not found (skipping): $src"
    fi
}

# 1. Move Big-BWT Logs (Original and Reverse)
# These are typically created at ${base}.log and ${base}.rev.log
safe_move_log "${base}.log" "$log_run_dir"
safe_move_log "${base}.rev.log" "$log_run_dir"

# 2. Move Moni-Align Logs
# Moni usually generates logs based on the output prefix if applicable
if [ -n "$moni_output" ]; then
    # Attempt to move the standard log file if moni produces one ending in .log
    safe_move_log "${moni_output}.log" "$log_run_dir"
    
    # If moni produces other log-like files starting with the prefix, try moving those too
    # Using a loop with nullglob to avoid errors if no files match
    shopt -s nullglob
    for logfile in "${moni_output}"*.log; do
        # Check to ensure we aren't moving the file we just moved
        if [ -f "$logfile" ]; then
             mv "$logfile" "$log_run_dir/"
             echo "Moved moni log: $logfile -> $log_run_dir"
        fi
    done
    shopt -u nullglob
fi

echo "Log reorganization complete."
echo "-------------------------------------------------------------"


# Cleanup
echo "Remove temporary files..."

# NOTE: "${base}.log" and "${base}.rev.log" were moved above, so they are safe.
rm -f "${base}.bwt" "${base}.rev.bwt" "${base}.ssa" "${base}.rev.ssa" \
      "${base}.esa" "${base}.rev.esa" \
      "${base}" "${base}.rev"
      
echo "Temporary files removed!" 

end_time=$(date +%s)
echo "Total time elapsed: $((end_time - start_time)) seconds."
echo "Script log saved to: $log_filename"