#!/bin/bash

#####################################
### Copy raw Cell Ranger matrices ###
#####################################

set -euo pipefail

# local folder containing Cell Ranger count output directories
src_base=~/cellranger/count
target_base=~/Thema_R/data/raw_counts

mkdir -p "$target_base"

echo "=== Copying raw feature matrices ==="

for sample_dir in "$src_base"/*; do
    [ -d "$sample_dir" ] || continue

    sample=$(basename "$sample_dir")

    raw_dir="$sample_dir/outs/raw_feature_bc_matrix"
    raw_h5="$sample_dir/outs/raw_feature_bc_matrix.h5"

    # skip non-sample folders
    if [ ! -d "$raw_dir" ] && [ ! -f "$raw_h5" ]; then
        echo "[SKIP] $sample is not a valid sample folder"
        continue
    fi

    # convert sample name:
    # 9-9NM   -> 9NM
    # 10-10TA -> 10TA
    left="${sample%%-*}"
    right="${sample#*-}"
    short_name="${left}${right#${left}}"

    # copy raw_feature_bc_matrix directory
    if [ -d "$raw_dir" ]; then
        out_dir="$target_base/${short_name}_raw_feature_bc_matrix"
        mkdir -p "$out_dir"

        cp "$raw_dir/barcodes.tsv.gz" "$out_dir/"
        cp "$raw_dir/features.tsv.gz" "$out_dir/"
        cp "$raw_dir/matrix.mtx.gz" "$out_dir/"

        echo "[OK] Copied folder matrix for $sample -> ${short_name}_raw_feature_bc_matrix"
    else
        echo "[WARN] Missing raw_feature_bc_matrix directory for $sample"
    fi

    # copy raw_feature_bc_matrix.h5
    if [ -f "$raw_h5" ]; then
        cp "$raw_h5" "$target_base/${short_name}_raw_feature_bc_matrix.h5"
        echo "[OK] Copied H5 matrix for $sample -> ${short_name}_raw_feature_bc_matrix.h5"
    else
        echo "[WARN] Missing raw_feature_bc_matrix.h5 for $sample"
    fi
done

echo "=== Completed successfully ==="
