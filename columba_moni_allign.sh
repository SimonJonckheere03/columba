#!/bin/bash

# ==============================================================================
# SCRIPT: Sequential Alignment Wrapper for Columba and Moni
# DESCRIPTION: Runs both aligners with user-definable parameters.
#              Organizes outputs into timestamped folders (Logs vs Results).
#
# USAGE: ./run_alignment_pipeline.sh -L <output_dir> -c <columba_idx> -m <moni_idx> ...
# ==============================================================================

# --- 1. DEFAULT VARIABLES -----------------------------------------------------

# Executables (Change these if they are not in the current folder)
COLUMBA_EXEC="./columba"
MONI_EXEC="./moni"

# Default Parameters (Can be overridden by flags)
COLUMBA_MODE="all"       # Default for Columba -a
COLUMBA_DIST="0"         # Default for Columba -e
THREADS=1                # Default threads

# Initialize paths as empty to check for missing arguments later
COLUMBA_INDEX=""
MONI_INDEX=""
READS_FILE=""
BASE_OUTPUT_DIR=""

# --- 2. USAGE FUNCTION --------------------------------------------------------

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Required:"
    echo "  -L <path>   Base directory for logs and results (Required)"
    echo "  -c <path>   Path to Columba Index"
    echo "  -m <path>   Path to Moni Index"
    echo "  -f <path>   Path to Reads file (.fastq)"
    echo ""
    echo "Optional:"
    echo "  -a <str>    Columba Alignment Mode (default: 'all')"
    echo "  -e <int>    Columba Max Distance (default: 0)"
    echo "  -t <int>    Number of threads (default: 1)"
    echo "  -h          Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -L ./runs -c ./idx/columba -m ./idx/moni -f ./data/reads.fq -a best"
    exit 1
}

# --- 3. PARSE COMMAND LINE ARGUMENTS ------------------------------------------

while getopts "L:c:m:f:a:e:t:h" opt; do
  case $opt in
    L) BASE_OUTPUT_DIR="$OPTARG" ;;
    c) COLUMBA_INDEX="$OPTARG" ;;
    m) MONI_INDEX="$OPTARG" ;;
    f) READS_FILE="$OPTARG" ;;
    a) COLUMBA_MODE="$OPTARG" ;;
    e) COLUMBA_DIST="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

# --- 4. VALIDATION ------------------------------------------------------------

# Check if required arguments are set
if [[ -z "$BASE_OUTPUT_DIR" || -z "$COLUMBA_INDEX" || -z "$MONI_INDEX" || -z "$READS_FILE" ]]; then
    echo "Error: Missing required arguments."
    echo "You must specify output dir (-L), indices (-c, -m), and reads (-f)."
    echo "------------------------------------------------------------"
    usage
fi

# --- 5. DIRECTORY SETUP -------------------------------------------------------

# Generate Timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Define run-specific folder structure
RUN_DIR="${BASE_OUTPUT_DIR}/run_${TIMESTAMP}"
LOGS_DIR="${RUN_DIR}/logs"
RESULTS_DIR="${RUN_DIR}/results"

echo "------------------------------------------------------------"
echo "Initializing Run: ${TIMESTAMP}"
echo "------------------------------------------------------------"
echo "Configuration:"
echo "  [Reads]        ${READS_FILE}"
echo "  [Columba Idx]  ${COLUMBA_INDEX}"
echo "  [Columba Mode] ${COLUMBA_MODE} (Dist: ${COLUMBA_DIST})"
echo "  [Moni Idx]     ${MONI_INDEX}"
echo "  [Output Dir]   ${RUN_DIR}"
echo "------------------------------------------------------------"

# Create folders
mkdir -p "${LOGS_DIR}"
mkdir -p "${RESULTS_DIR}"

# Define Output Filenames
COLUMBA_LOG="${LOGS_DIR}/columba.log"
COLUMBA_RESULT="${RESULTS_DIR}/columba_output.sam"

MONI_LOG="${LOGS_DIR}/moni.log"
MONI_RESULT="${RESULTS_DIR}/moni_output.sam"

# --- 6. EXECUTION: COLUMBA ----------------------------------------------------

echo "Step 1/2: Running Columba..."
echo "  > Saving log to: ${COLUMBA_LOG}"

# Define a temporary file for stats inside the logs folder
MEM_STATS="${LOGS_DIR}/columba_resource_stats.tmp"

# Run Columba wrapped in /usr/bin/time.
# We remove the '-o' flag and instead use '2>' to force-redirect 
# all stderr (which includes the timing stats) to our temp file.
/usr/bin/time -v \
    $COLUMBA_EXEC \
    -r "${COLUMBA_INDEX}" \
    -f "${READS_FILE}" \
    -a "${COLUMBA_MODE}" \
    -e "${COLUMBA_DIST}" \
    -t "${THREADS}" \
    -l "${COLUMBA_LOG}" \
    -o "${COLUMBA_RESULT}" \
    2> "${MEM_STATS}"

STATUS=$?

# Append the captured stats to the main log file
if [ -s "${MEM_STATS}" ]; then
    echo "" >> "${COLUMBA_LOG}"
    echo "==============================" >> "${COLUMBA_LOG}"
    echo "   SYSTEM RESOURCE USAGE      " >> "${COLUMBA_LOG}"
    echo "==============================" >> "${COLUMBA_LOG}"
    cat "${MEM_STATS}" >> "${COLUMBA_LOG}"
    rm "${MEM_STATS}"
else
    echo "Warning: No resource usage stats were captured." >> "${COLUMBA_LOG}"
fi

if [ $STATUS -eq 0 ]; then
    echo "  > Columba finished successfully."
else
    echo "  > ERROR: Columba failed. Check log for details."
    # exit 1 
fi

# --- 7. EXECUTION: MONI -------------------------------------------------------

echo "------------------------------------------------------------"
echo "Step 2/2: Running Moni..."
echo "  > Saving log to: ${MONI_LOG}"

# Moni writes logs to stdout/stderr, so we pipe everything to the log file.
# To add timing, we wrap the whole command.
/usr/bin/time -v \
    $MONI_EXEC align \
    -i "${MONI_INDEX}" \
    -p "${READS_FILE}" \
    -o "${MONI_RESULT}" \
    > "${MONI_LOG}" 2>&1

# Note: Since we redirected Moni's output to the file using '>', 
# the 'time' output (which goes to stderr) will also end up in the 
# log file automatically because of the '2>&1' at the end.

if [ $? -eq 0 ]; then
    echo "  > Moni finished successfully."
else
    echo "  > ERROR: Moni failed. Check log for details."
fi

# ==============================================================================

echo "------------------------------------------------------------"
echo "Pipeline complete."
echo "Results stored in: ${RUN_DIR}"