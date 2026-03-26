#!/bin/bash

# Script: columba_build_pfp.sh
# Description: Columba build process using the Moni Hybrid Pipeline with VCF support.
# Author: Lore Depuydt - lore.depuydt@ugent.be
# Modified by Simon Jonckheere - simon.jonckheere@ugent.be

# Capture start time
start_time=$(date +%s)

columba_build_exe="./columba_build"

# --- EXECUTABLE PLACEHOLDERS ---
pfp_exe="./bin/pfp++"
pfp64_exe="./bin/pfp++64"
bwtparse_exe="../external/Big-BWT/bwtparse"
bwtparse64_exe="../external/Big-BWT/bwtparse64"
pfbwtNT_exe=".../external/Big-BWT/pfbwtNT.x"
pfbwtNT64_exe="../external/Big-BWT/pfbwtNT64.x"
# -------------------------------

seedLength=100
ws=0
mod=0
haplotype="1"
vcf_file=""
fasta_files=()

showUsage() {
    echo "Usage: $0 [-l <seedLength>] [-w <ws>] [-p <mod>] -r <index_name> [-f <fasta_file.fa.gz>] [-V <vcf_file>] [-H <haplotype>]"
    echo
    echo "Required arguments:"
    echo "  -r <index_name>       Name/location of the index to be created."
    echo
    echo "Optional arguments:"
    echo "  -f <fasta_file>       The reference FASTA file (already in .fa.gz format)."
    echo "  -F <fasta_file_list>  Path to a file containing a list of FASTA files (ensure pfp++ supports multiple if used)."
    echo "  -V <vcf_file>         VCF file (bgzipped and indexed) to apply to the reference."
    echo "  -H <haplotype>        Haplotype to extract from VCF (default: 1)."
    echo "  -l <seedLength>       Seed length for replacing non-ACGT characters (default: $seedLength). *May be obsolete if preprocess is removed.*"
    echo "  -w <ws>               Window size for Big-BWT."
    echo "  -p <mod>              Mod value for Big-BWT."
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


# Parse command-line options
while getopts ":l:r:f:F:w:p:V:H:" opt; do
    case $opt in
        l) seedLength=$OPTARG ;;
        r) index_name=$OPTARG ;;
        f)
            fasta_files+=("$OPTARG")
            while [[ $OPTIND -le $# && ! ${!OPTIND} =~ ^- ]]; do
                fasta_files+=("${!OPTIND}")
                OPTIND=$((OPTIND + 1))
            done
            ;;
        F)
            if [[ -f $OPTARG ]]; then
                while IFS= read -r line; do
                    if [[ -n "$line" ]]; then
                        fasta_files+=("$line")
                    fi
                done <"$OPTARG"
            else
                echo "Error: File '$OPTARG' not found." >&2
                exit 1
            fi
            ;;
        V) vcf_file=$OPTARG ;;
        H) haplotype=$OPTARG ;;
        w) ws=$OPTARG ;;
        p) mod=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; showUsage; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; showUsage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ -z "$index_name" ] || [ "${#fasta_files[@]}" -eq 0 ]; then
    showUsage
    exit 1
fi

echo "Welcome to the Columba build process with bidirectional prefix-free parsing!"
echo "-------------------------------------------------------------"

base="${index_name}"
# Assuming the user passes a single .fa.gz file as the reference
ref_archive="${fasta_files[0]}"

# 1. Run PFP++ ONCE (Generates both Forward and Reverse parses)
echo "Start Bidirectional Prefix-Free Parsing via pfp++..."
w_val="${ws:-10}"
p_val="${mod:-100}"
if [ "$w_val" -eq 0 ]; then w_val=10; fi
if [ "$p_val" -eq 0 ]; then p_val=100; fi

# Passing the input .fa.gz directly to -r
runCommandWithTime "$pfp64_exe" -w "$w_val" -p "$p_val" -o "$base" -c -i -l --acgt-only -r "$ref_archive" -v "$vcf_file" -H "$haplotype"
echo "Prefix-free parsing done! (Generated .parse and .rev.parse)"
echo "-------------------------------------------------------------"

# 2. Run BigBWT on the FORWARD strings
echo "Running Big-BWT on FORWARD parse..."
runCommandWithTime "$bwtparse64_exe" "$base" -s
runCommandWithTime "$pfbwtNT64_exe" -w "$w_val" -s -e "$base"

# 3. Run BigBWT on the REVERSE strings
echo "Running Big-BWT on REVERSE parse..."
runCommandWithTime "$bwtparse64_exe" "${base}.rev" -s
runCommandWithTime "$pfbwtNT64_exe" -w "$w_val" -s -e "${base}.rev"
echo "Big-BWT generation complete."
echo "-------------------------------------------------------------"

# 4. Build Columba Index
echo "Start building the Columba index..."
runCommandWithTime "$columba_build_exe" --pfp -r "$index_name"
echo "Columba index built!"
echo "-------------------------------------------------------------"

# 5. Cleanup
echo "Cleaning up temporary files..."
rm -f "${base}.bwt" "${base}.rev.bwt" "${base}.ssa" "${base}.rev.ssa" "${base}.esa" "${base}.rev.esa"
rm -f "${base}.log" "${base}.rev.log"
rm -f "${base}.parse" "${base}.dict" "${base}.dicz" "${base}.occ" "${base}.ilist" "${base}.last" "${base}.bwlast"
rm -f "${base}.rev.parse" "${base}.rev.dict" "${base}.rev.dicz" "${base}.rev.occ" "${base}.rev.ilist" "${base}.rev.last" "${base}.rev.bwlast"

echo "Temporary files removed!"
echo "Total time elapsed: $(($(date +%s) - start_time)) seconds."