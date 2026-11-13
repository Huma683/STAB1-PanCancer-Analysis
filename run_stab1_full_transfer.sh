#!/bin/bash
set -e

echo "===== STAB1 Condor transfer job started ====="
hostname
date
echo "Condor scratch directory:"
pwd
RSCRIPT=$(which Rscript)
echo "Using Rscript:"
echo $RSCRIPT
$RSCRIPT --version

echo "Files transferred into scratch:"
ls -lh

mkdir -p data scripts results logs

echo "Preparing input files..."
cp "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz" data/
cp "Survival_SupplementalTable_S1_20171025_xena_sp" data/
cp "stab1_baseR_full_survival_stage_condor.R" scripts/

echo "Running STAB1 full analysis in Condor scratch..."
$RSCRIPT scripts/stab1_baseR_full_survival_stage_condor.R

echo "Checking generated outputs:"
find results -type f -maxdepth 3 -print
find results -type f -maxdepth 3 -exec ls -lh {} \;

echo "Compressing all results..."
tar -czf stab1_full_condor_results.tar.gz results

echo "Final output file:"
ls -lh stab1_full_condor_results.tar.gz

echo "===== STAB1 Condor transfer job finished successfully ====="
date
