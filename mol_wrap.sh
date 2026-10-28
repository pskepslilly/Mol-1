#!/bin/bash

print_help() {
        cat << 'EOF'
Usage:
    mol_wrap.sh -i <input_fasta> -o <output_dir> [-s <num_seeds>] [mol options]

Wrapper options (intercepted here):
    -i, --input (str):       [Required] Input FASTA file.
    -o, --out_dir (str):     [Default: ./output] Output directory for results.
    -s, --seeds (int):       [Default: 5] Number of seeds to generate for the model.
    -b, --base_seed (int):   [Default: 2026] The base seed for seeds to be built off of
    -t, --template:          Whether to run template search and include templates in prediction.

Forwarded mol options:
    -c, --cycle (int):       [Default: 10] Number of Pairformer cycles.
    -p, --step (int):        [Default: 200] Number of diffusion steps.
    -e, --sample (int):      [Default: 5] Samples per seed.
    -d, --dtype (str):       [Default: bf16] Inference data type (bf16, fp32).
    -n, --model_name (str):  [Default: mol-1-base] Model name.
    --use_msa (bool):        Whether to use protein MSA features.
    --use_default_params:    Auto-load recommended defaults for the model.
    --trimul_kernel (str):   Triangle multiplicative kernel ('cuequivariance', 'torch').
    --triatt_kernel (str):   Triangle attention kernel ('triattention', 'cuequivariance', etc.).
    --use_rna_msa (bool):    Enable RNA MSA features (v1.0.0+ only).
EOF
}

fasta_to_mol_json() {
    local fasta_file="$1"
    local output_json="${2:-output.json}"
    local msa_base_dir="${3:-./msas}"
    local name="${4:-$(basename "$fasta_file" .fasta)}"
    local num_seeds="${5:-5}"
    local base_seed="${6:-2026}"

    if [ ! -f "$fasta_file" ]; then
        echo "Error: FASTA file not found: $fasta_file"
        return 1
    fi

    python3 - "$fasta_file" "$output_json" "$msa_base_dir" "$name" "$num_seeds" "$base_seed" << 'PYTHON_EOF'
import json
import sys
import re
from collections import OrderedDict

fasta_file = sys.argv[1]
output_json = sys.argv[2]
msa_base_dir = sys.argv[3]
name = sys.argv[4]
num_seeds = int(sys.argv[5])
BASE_SEED = int(sys.argv[6])

def extract_count_from_header(header):
    """Extract explicitly specified count from header"""
    patterns = [
        r'count[=:](\d+)', r'copies[=:](\d+)', r'stoich[=:](\d+)',
        r'_count_(\d+)', r'\(count=(\d+)\)', r'\[(\d+)x\]',
        r'\sx(\d+)$', r'_x(\d+)$',
    ]
    for pattern in patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None  # Return None if not explicitly specified

# Parse FASTA file
all_sequences = []
current_seq = ""
current_header = ""

with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                all_sequences.append({
                    'header': current_header,
                    'sequence': current_seq
                })
            current_seq = ""
            current_header = line
        else:
            current_seq += line.replace(" ", "").replace("\n", "")

# Don't forget last sequence
if current_seq:
    all_sequences.append({
        'header': current_header,
        'sequence': current_seq
    })

# Group identical sequences and count them
# Use OrderedDict to preserve order of first occurrence
unique_sequences = OrderedDict()

for entry in all_sequences:
    seq = entry['sequence']
    header = entry['header']

    # Check if this sequence already exists
    if seq in unique_sequences:
        # Check if count was explicitly specified in header
        explicit_count = extract_count_from_header(header)
        if explicit_count is not None:
            # Use explicit count, don't auto-increment
            # This handles cases like ">chain1|count=2" appearing twice
            unique_sequences[seq]['explicit_counts'].append(explicit_count)
        else:
            # Auto-increment count for identical sequences
            unique_sequences[seq]['auto_count'] += 1
    else:
        # First occurrence of this sequence
        explicit_count = extract_count_from_header(header)
        unique_sequences[seq] = {
            'sequence': seq,
            'first_header': header,
            'auto_count': 1,
            'explicit_counts': [explicit_count] if explicit_count is not None else []
        }

# Build final sequences list
sequences = []
seq_num = 0

for seq, info in unique_sequences.items():

    # Determine final count
    if info['explicit_counts']:
        # If any explicit counts were specified, use the sum
        count = sum(info['explicit_counts'])
        count_source = f"explicit: {'+'.join(map(str, info['explicit_counts']))}"
    else:
        # Use auto-detected count
        count = info['auto_count']
        count_source = f"auto-detected ({info['auto_count']} identical sequences)"

    sequences.append({
        "proteinChain": {
            "sequence": seq,
            "count": count,
            "pairedMsaPath":f"{msa_base_dir}/{seq_num}/pairing.a3m",
            "unpairedMsaPath":f"{msa_base_dir}/{seq_num}/non_pairing.a3m"
        }
    })

    print(f"Chain {seq_num}: count={count} ({count_source})", file=sys.stderr)
    seq_num += 1

# Generate seeds
model_seeds = [BASE_SEED + i for i in range(num_seeds)]

output = [{
    "sequences": sequences,
    "modelSeeds": model_seeds,
    "name": name
}]

with open(output_json, "w") as f:
    json.dump(output, f, indent=4)

print(f"\n✓ Created {output_json}:", file=sys.stderr)
print(f"  Unique chains: {len(sequences)}", file=sys.stderr)
total_copies = sum(seq["proteinChain"]["count"] for seq in sequences)
print(f"  Total copies: {total_copies}", file=sys.stderr)
print(f"  Seeds: [{model_seeds[0]}..{model_seeds[-1]}] (n={len(model_seeds)})", file=sys.stderr)
PYTHON_EOF
}

concatenate_fasta_for_multimer() {
    local input_fasta="$1"
    local output_fasta="${2:-${input_fasta%.fasta}_multimer.fasta}"

    if [ ! -f "$input_fasta" ]; then
        echo "Error: Input FASTA not found: $input_fasta" >&2
        return 1
    fi

    # Extract all sequences and join with colons
    local sequences=()
    local current_seq=""

    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            if [ -n "$current_seq" ]; then
                sequences+=("$current_seq")
            fi
            current_seq=""
        else
            # Remove whitespace and append
            current_seq="${current_seq}${line//[[:space:]]/}"
        fi
    done < "$input_fasta"

    # Don't forget last sequence
    if [ -n "$current_seq" ]; then
        sequences+=("$current_seq")
    fi

    # Join with colons
    local concatenated=$(IFS=:; echo "${sequences[*]}")

    # Write new FASTA
    local header=">multimer_${#sequences[@]}_chains"
    echo "$header" > "$output_fasta"
    echo "$concatenated" >> "$output_fasta"

    echo "Created multimer FASTA: $output_fasta (${#sequences[@]} chains)" >&2
    echo "$output_fasta"
}

BUILD_PATH="/trainfs/data/comp_bio/apps/Mol-1"
VENV="${BUILD_PATH}/venv"
PROF="/trainfs/data/comp_bio/conda.profile"

source ${PROF}
conda activate ${VENV}

export MOL_ROOT_DIR="/trainfs/data/comp_bio/apps/Mol-1/data"

COLAB_DB="/trainfs/data/comp_bio/data/COLAB_DB"
COLAB_SEARCH="/trainfs/data/comp_bio/bin/colabfold_search"

if [[ "$1" == "-h" || "$1" == "--help" || "$#" -eq 0 ]]; then
    print_help
    exit 0
fi

INPUT_FASTA=""
OUTPUT_DIR=""
NUM_SEEDS="5"
USE_TEMPLATES="0"
BASE_SEED="2026"

EXTRA_ARGS=()
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h|--help)
            print_help
            exit 0
            ;;
        -i|--input)
            INPUT_FASTA="$2"
            shift 2
            ;;
        -o|--out_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -b|--base_seed)
            BASE_SEED="$2"
            shift 2
            ;;
        -s|--seeds)
            NUM_SEEDS="$2"
            shift 2
            ;;
        -t|--template)
            echo "Using templates during the prediction"
            USE_TEMPLATES="1"
            EXTRA_ARGS+=("--use_template" "true")
            shift 1
            ;;
        --use_seeds_in_json)
            echo "Skipping disallowed option: $1" >&2
            if [[ -n "$2" && "$2" != -* ]]; then
                shift 2
            else
                shift 1
            fi
            ;;
        --*)
            EXTRA_ARGS+=("$1")
            shift 1
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift 1
            ;;
    esac
done

if [[ -z "$INPUT_FASTA" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: input FASTA and output directory are required." >&2
    print_help
    exit 1
fi

RUN_NAME=${INPUT_FASTA##*/}
RUN_NAME=${RUN_NAME%.*}

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

MULTIMER_FASTA="${OUTPUT_DIR}/${RUN_NAME}_multimer.fasta"
concatenate_fasta_for_multimer "$INPUT_FASTA" "$MULTIMER_FASTA"

python3 ${BUILD_PATH}/scripts/colabfold_msa.py ${MULTIMER_FASTA} ${COLAB_DB} ${OUTPUT_DIR} --colabsearch ${COLAB_SEARCH} --db1 uniref30_2302_db --db2 pdb100_230517 --db3 colabfold_envdb_202108_db --gpu 1

fasta_to_mol_json ${INPUT_FASTA} ${OUTPUT_DIR}/input.json ${OUTPUT_DIR}/msa "${RUN_NAME}" ${NUM_SEEDS} ${BASE_SEED}

INPUT_JSON="${OUTPUT_DIR}/input.json"

mol pred \
    -i ${INPUT_JSON} \
    -o ${OUTPUT_DIR}/predictions \
    --use_seeds_in_json true \
    "${EXTRA_ARGS[@]}" \
