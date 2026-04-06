#!/bin/bash
set -euo pipefail

# Script: columba_build_pfp.sh
# Description: Columba build process using the Moni Hybrid Pipeline with VCF support.
# Author: Lore Depuydt - lore.depuydt@ugent.be
# Modified by Simon Jonckheere - simon.jonckheere@ugent.be

start_time=$(date +%s)
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -x "${script_dir}/bin/columba_build" ]]; then
    build_dir="${script_dir}"
elif [[ -x "${script_dir}/../../build_RLC/bin/columba_build" ]]; then
    build_dir=$(cd "${script_dir}/../../build_RLC" && pwd)
else
    echo "Error: Could not locate the Columba build directory from $script_dir" >&2
    exit 1
fi
repo_dir=$(cd "${build_dir}/.." && pwd)

columba_build_exe="${build_dir}/bin/columba_build"
pfp64_exe="${build_dir}/bin/pfp++64"
bwt_newscan_exe="${repo_dir}/external/Big-BWT/newscanNT.x"
bwtparse64_exe="${repo_dir}/external/Big-BWT/bwtparse64"
pfbwtNT64_exe="${repo_dir}/external/Big-BWT/pfbwtNT64.x"

seedLength=100
ws=0
mod=0
haplotype="1"
max_samples=""
samples_file=""
name_prefix=""
vcf_file=""
fasta_files=()
keep_intermediates=0

showUsage() {
    echo "Usage: $0 [-l <seedLength>] [-w <ws>] [-p <mod>] -r <index_name> [-f <fasta_file.fa.gz>] [-V <vcf_file>] [-H <haplotype>] [-M <max_samples>] [-S <samples_file>] [-P <name_prefix>] [-k]"
    echo
    echo "Required arguments:"
    echo "  -r <index_name>       Name/location of the index to be created."
    echo
    echo "Optional arguments:"
    echo "  -f <fasta_file>       The reference FASTA file (already in .fa.gz format)."
    echo "  -F <fasta_file_list>  Path to a file containing a list of FASTA files."
    echo "  -V <vcf_file>         VCF file (bgzipped and indexed) to apply to the reference."
    echo "  -H <haplotype>        Haplotype to extract from VCF (default: $haplotype)."
    echo "  -M <max_samples>      Maximum number of VCF samples to parse."
    echo "  -S <samples_file>     File containing VCF sample IDs to parse."
    echo "  -P <name_prefix>      Prefix added to generated haplotype sequence names."
    echo "  -l <seedLength>       Seed length for Columba-compatible non-ACGT replacement (default: $seedLength)."
    echo "  -w <ws>               Window size for PFP / Big-BWT."
    echo "  -p <mod>              Mod value for PFP / Big-BWT."
    echo "  -k                    Keep PFP and Big-BWT intermediate artifacts after success."
}

runCommandWithTime() {
    local command="$1"
    shift
    (/usr/bin/time -v "$command" "$@") || {
        local status=$?
        echo "Error: Command '$command $*' failed with exit status $status." >&2
        exit $status
    }
}

requireExecutable() {
    local path="$1"
    if [[ ! -x "$path" ]]; then
        echo "Error: Required executable not found or not executable: $path" >&2
        exit 1
    fi
}

requireFile() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "Error: Required file not found: $path" >&2
        exit 1
    fi
}

requireNonEmptyFile() {
    local path="$1"
    requireFile "$path"
    if [[ ! -s "$path" ]]; then
        echo "Error: Required file is empty: $path" >&2
        exit 1
    fi
}

requireArtifacts() {
    local label="$1"
    shift
    echo "Validating $label artifacts..."
    local artifact
    for artifact in "$@"; do
        requireNonEmptyFile "$artifact"
    done
}

writeSequenceMetadataFromPFP() {
    local base="$1"
    local ref_archive="$2"
    python3 - "$base" "$ref_archive" <<'PY'
import gzip
import struct
import sys
from pathlib import Path

base = Path(sys.argv[1])
ref_archive = Path(sys.argv[2])
lidx_path = Path(f"{base}.lidx")

if not lidx_path.is_file():
    raise SystemExit(f"Missing .lidx file: {lidx_path}")
if not ref_archive.is_file():
    raise SystemExit(f"Missing reference FASTA: {ref_archive}")

names = []
payload_lengths = []
for raw_line in lidx_path.read_text().splitlines():
    line = raw_line.strip()
    if not line:
        continue
    name, length = line.rsplit(maxsplit=1)
    names.append(name)
    payload_lengths.append(int(length))

positions = []
payload_pos = 0
for payload_len in payload_lengths:
    positions.append(payload_pos)
    payload_pos += payload_len
positions.append(payload_pos)

(base.parent / f"{base.name}.pos").write_bytes(
    b"".join(struct.pack("<Q", value) for value in positions)
)

with (base.parent / f"{base.name}.sna").open("wb") as handle:
    for name in names:
        encoded = name.encode()
        handle.write(struct.pack("<Q", len(encoded)))
        handle.write(encoded)

(base.parent / f"{base.name}.fsid").write_bytes(struct.pack("<Q", 0))

reference_names = []
reference_lengths = []
with gzip.open(ref_archive, "rt") as handle:
    current_name = None
    current_length = 0
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                reference_names.append(current_name)
                reference_lengths.append(current_length)
            current_name = line[1:].split()[0]
            current_length = 0
        else:
            current_length += len(line)
    if current_name is not None:
        reference_names.append(current_name)
        reference_lengths.append(current_length)

with (base.parent / f"{base.name}.headerSN.bin").open("wb") as handle:
    for name, length in zip(reference_names, reference_lengths):
        handle.write(f"@SQ\tSN:{name}\tLN:{length}\n".encode())
PY
}

cleanupArtifacts() {
    local base="$1"
    rm -f "${base}.bwt" "${base}.rev.bwt" "${base}.ssa" "${base}.rev.ssa" "${base}.esa" "${base}.rev.esa"
    rm -f "${base}.log" "${base}.rev.log"
    rm -f "${base}.parse" "${base}.dict" "${base}.dicz" "${base}.occ" "${base}.ilist" "${base}.last" "${base}.bwlast" "${base}.sai" "${base}.bwsai"
    rm -f "${base}.rev.parse" "${base}.rev.dict" "${base}.rev.dicz" "${base}.rev.occ" "${base}.rev.ilist" "${base}.rev.last" "${base}.rev.bwlast" "${base}.rev.sai" "${base}.rev.bwsai"
    rm -f "${base}.parse.u64" "${base}.occ.u64" "${base}.rev.parse.u64" "${base}.rev.occ.u64"
    rm -f "${base}.bbwt" "${base}.bbwt.map" "${base}.bbwt.parse" "${base}.bbwt.dict" "${base}.bbwt.occ" "${base}.bbwt.last" "${base}.bbwt.sai" "${base}.bbwt.ilist" "${base}.bbwt.bwlast" "${base}.bbwt.bwsai" "${base}.bbwt.bwt" "${base}.bbwt.ssa" "${base}.bbwt.esa"
    rm -f "${base}.rev.bbwt" "${base}.rev.bbwt.map" "${base}.rev.bbwt.parse" "${base}.rev.bbwt.dict" "${base}.rev.bbwt.occ" "${base}.rev.bbwt.last" "${base}.rev.bbwt.sai" "${base}.rev.bbwt.ilist" "${base}.rev.bbwt.bwlast" "${base}.rev.bbwt.bwsai" "${base}.rev.bbwt.bwt" "${base}.rev.bbwt.ssa" "${base}.rev.bbwt.esa"
}

installBigBWTOutputs() {
    local bbwt_prefix="$1"
    local target_prefix="$2"
    mv "${bbwt_prefix}.bwt" "${target_prefix}.bwt"
    mv "${bbwt_prefix}.ssa" "${target_prefix}.ssa"
    mv "${bbwt_prefix}.esa" "${target_prefix}.esa"
}

while getopts ":l:r:f:F:w:p:V:H:M:S:P:k" opt; do
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
        M) max_samples=$OPTARG ;;
        S) samples_file=$OPTARG ;;
        P) name_prefix=$OPTARG ;;
        w) ws=$OPTARG ;;
        p) mod=$OPTARG ;;
        k) keep_intermediates=1 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; showUsage; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; showUsage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [[ -z "${index_name:-}" || ${#fasta_files[@]} -eq 0 ]]; then
    showUsage
    exit 1
fi

if [[ ${#fasta_files[@]} -ne 1 ]]; then
    echo "Error: This script currently supports exactly one FASTA archive." >&2
    exit 1
fi

base="${index_name}"
ref_archive="${fasta_files[0]}"
mkdir -p "$(dirname "$base")"

requireExecutable "$columba_build_exe"
requireExecutable "$pfp64_exe"
requireExecutable "$bwt_newscan_exe"
requireExecutable "$bwtparse64_exe"
requireExecutable "$pfbwtNT64_exe"
requireFile "$ref_archive"
if [[ -n "$vcf_file" ]]; then
    requireFile "$vcf_file"
fi
if [[ -n "$samples_file" ]]; then
    requireFile "$samples_file"
fi
if [[ -n "$max_samples" ]]; then
    if ! [[ "$max_samples" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: -M <max_samples> must be a positive integer." >&2
        exit 1
    fi
fi
if ! [[ "$seedLength" =~ ^[0-9]+$ ]]; then
    echo "Error: -l <seedLength> must be a non-negative integer." >&2
    exit 1
fi

w_val="${ws:-10}"
p_val="${mod:-100}"
if [[ "$w_val" -eq 0 ]]; then w_val=10; fi
if [[ "$p_val" -eq 0 ]]; then p_val=100; fi
if [[ "$w_val" -lt 4 ]]; then
    echo "Error: Big-BWT requires a window size of at least 4." >&2
    exit 1
fi

echo "Welcome to the Columba build process with bidirectional prefix-free parsing!"
echo "-------------------------------------------------------------"
echo "Build directory: $build_dir"
echo "Repository directory: $repo_dir"
echo "Output base: $base"

echo "Start Bidirectional Prefix-Free Parsing via pfp++..."
pfp_cmd=("$pfp64_exe" -w "$w_val" -p "$p_val" -o "$base" -c -i -l --acgt-only --seed-length "$seedLength" -r "$ref_archive")
if [[ -n "$vcf_file" ]]; then
    pfp_cmd+=(-v "$vcf_file" -H "$haplotype")
    if [[ -n "$max_samples" ]]; then
        pfp_cmd+=(-m "$max_samples")
    fi
    if [[ -n "$samples_file" ]]; then
        pfp_cmd+=(-S "$samples_file")
    fi
    if [[ -n "$name_prefix" ]]; then
        pfp_cmd+=(-P "$name_prefix")
    fi
fi
runCommandWithTime "${pfp_cmd[@]}"
requireArtifacts "forward PFP" \
    "${base}.parse" "${base}.dict" "${base}.last" "${base}.sai" "${base}.occ" "${base}.lidx" "${base}.ldx"
requireArtifacts "Big-BWT payloads" \
    "${base}.bbwt" "${base}.rev.bbwt"
writeSequenceMetadataFromPFP "$base" "$ref_archive"
requireArtifacts "sequence metadata" \
    "${base}.headerSN.bin" "${base}.pos" "${base}.sna" "${base}.fsid"
echo "Prefix-free parsing done."
echo "-------------------------------------------------------------"

parser_mod="$p_val"
if [[ "$parser_mod" -lt 10 ]]; then parser_mod=10; fi

echo "Running Big-BWT on FORWARD payload..."
runCommandWithTime "$bwt_newscan_exe" "${base}.bbwt" -w "$w_val" -p "$parser_mod" -s
requireArtifacts "forward payload parser" \
    "${base}.bbwt.parse" "${base}.bbwt.dict" "${base}.bbwt.last" "${base}.bbwt.sai" "${base}.bbwt.occ"
runCommandWithTime "$bwtparse64_exe" "${base}.bbwt" -s
requireArtifacts "forward payload bwtparse" \
    "${base}.bbwt.ilist" "${base}.bbwt.bwlast" "${base}.bbwt.bwsai"
runCommandWithTime "$pfbwtNT64_exe" -w "$w_val" -s -e "${base}.bbwt"
requireArtifacts "forward payload Big-BWT" \
    "${base}.bbwt.bwt" "${base}.bbwt.ssa" "${base}.bbwt.esa"
installBigBWTOutputs "${base}.bbwt" "$base"
requireArtifacts "forward installed Big-BWT" \
    "${base}.bwt" "${base}.ssa" "${base}.esa"

echo "Running Big-BWT on REVERSE payload..."
runCommandWithTime "$bwt_newscan_exe" "${base}.rev.bbwt" -w "$w_val" -p "$parser_mod" -s
requireArtifacts "reverse payload parser" \
    "${base}.rev.bbwt.parse" "${base}.rev.bbwt.dict" "${base}.rev.bbwt.last" "${base}.rev.bbwt.sai" "${base}.rev.bbwt.occ"
runCommandWithTime "$bwtparse64_exe" "${base}.rev.bbwt" -s
requireArtifacts "reverse payload bwtparse" \
    "${base}.rev.bbwt.ilist" "${base}.rev.bbwt.bwlast" "${base}.rev.bbwt.bwsai"
runCommandWithTime "$pfbwtNT64_exe" -w "$w_val" -s -e "${base}.rev.bbwt"
requireArtifacts "reverse payload Big-BWT" \
    "${base}.rev.bbwt.bwt" "${base}.rev.bbwt.ssa" "${base}.rev.bbwt.esa"
installBigBWTOutputs "${base}.rev.bbwt" "${base}.rev"
requireArtifacts "reverse installed Big-BWT" \
    "${base}.rev.bwt" "${base}.rev.ssa" "${base}.rev.esa"
echo "Big-BWT generation complete."
echo "-------------------------------------------------------------"

echo "Start building the Columba index..."
runCommandWithTime "$columba_build_exe" --pfp -r "$index_name"
echo "Columba index built!"
echo "-------------------------------------------------------------"

if [[ "$keep_intermediates" -eq 1 ]]; then
    echo "Keeping intermediate PFP and Big-BWT artifacts."
else
    echo "Cleaning up temporary files..."
    cleanupArtifacts "$base"
    echo "Temporary files removed!"
fi

echo "Total time elapsed: $(($(date +%s) - start_time)) seconds."
