#!/bin/bash
set -e

echo "===== STAB1 Condor job started ====="
hostname
date

cd "$(dirname "$0")"

RSCRIPT=$(which Rscript)
echo "Using Rscript:"
echo $RSCRIPT

if [ ! -x "$RSCRIPT" ]; then
    echo "ERROR: Rscript not found at $RSCRIPT"
    exit 127
fi

$RSCRIPT --version

echo "Testing write permission..."
mkdir -p results/full_condor
echo "condor write test" > results/full_condor/condor_write_test.txt
cat results/full_condor/condor_write_test.txt
rm results/full_condor/condor_write_test.txt

echo "Running STAB1 full survival and stage-adjusted analysis through Condor..."
$RSCRIPT scripts/stab1_baseR_full_survival_stage_condor.R

echo "===== STAB1 Condor job finished successfully ====="
date
