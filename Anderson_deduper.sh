#!/bin/bash
#SBATCH --job-name=deduper
#SBATCH --output=deduper_%j.out
#SBATCH --error=deduper_%j.err
#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=5:00:00

# ------------------------------------------
# Reference-Based PCR Duplicate Removal
# Script: Anderson_deduper.py
# Cluster: BGMP (Talapas)
# ------------------------------------------

# ---- File paths ----
IN_SAM=/projects/bgmp/ica/bioinfo/Bi624/Deduper-iananderson221/C1_SE_uniqAlign.sorted.sam      
OUT_SAM=/projects/bgmp/ica/bioinfo/Bi624/Deduper-iananderson221/deduper_output_sam     
UMI_FILE=/projects/bgmp/ica/bioinfo/Bi624/Deduper-iananderson221/STL96.txt          
SCRIPT=/projects/bgmp/ica/bioinfo/Bi624/Deduper-iananderson221/Anderson_deduper.py
  


echo "=========================================="
echo "Deduper Job Started"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "------------------------------------------"
echo "Script:     $SCRIPT"
echo "Input SAM:  $IN_SAM"
echo "UMI list:   $UMI_FILE"
echo "Output SAM: $OUT_SAM"
echo "CPUs:       $SLURM_CPUS_ON_NODE"
echo "Mem:        16G"
echo "Time limit: 5:00:00"
echo "=========================================="

echo
echo ">>> Printing Anderson_deduper.py help (-h) to job log:"
echo "------------------------------------------"
python3 "$SCRIPT" -h
echo "------------------------------------------"
echo

echo ">>> Running deduplication..."
/usr/bin/time -v python3 "$SCRIPT" -f "$IN_SAM" -o "$OUT_SAM" -u "$UMI_FILE"

echo "------------------------------------------"
echo "Job finished: $(date)"
echo "=========================================="

