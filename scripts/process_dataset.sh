#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path_to_osh_pbf> [base_year]"
    echo "Example: $0 dataset/bhutan-internal.osh.pbf 2017"
    exit 1
fi

PBF_FILE=$1
BASE_YEAR=${2:-2017}

# 自动解析提取名字，比如从 bhutan-internal.osh.pbf 中提取出 bhutan
FILENAME=$(basename "$PBF_FILE")
DATASET_NAME="${FILENAME%-internal.osh.pbf}"
DATASET_NAME="${DATASET_NAME%.osh.pbf}" # 兜底逻辑

# 设定输出路径
WORK_DIR=$(dirname "$PBF_FILE")
PARSED_DIR="$WORK_DIR/${DATASET_NAME}_parsed"
WORKLOAD_DIR="$WORK_DIR/${DATASET_NAME}_workload"

echo "========================================"
echo " 🗄️ Dataset: $DATASET_NAME"
echo " 📁 PBF File: $PBF_FILE"
echo " 📅 Base Year: $BASE_YEAR"
echo "========================================"

echo ""
echo "⏳ [1/2] Extracting history to yearly CSVs..."
python3 scripts/extract_osm_history.py --input "$PBF_FILE" --outdir "$PARSED_DIR"

echo ""
echo "⚙️ [2/2] Running state machine and preparing workload..."
python3 scripts/prepare_workload.py --indir "$PARSED_DIR" --outdir "$WORKLOAD_DIR" --base_year "$BASE_YEAR"

echo ""
echo "========================================"
echo "✅ Success! Workload prepared at:"
echo "   $WORKLOAD_DIR"
echo "========================================"
echo "🎯 You can now run the benchmark:"
echo "./verify_bench -algo all -q_step 1000 -dir $WORKLOAD_DIR -start_year $((BASE_YEAR+1)) -end_year 2026"
