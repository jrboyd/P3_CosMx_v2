#!/bin/bash
SCRIPT=scripts/03_run_wrapper.sh 
echo $SCRIPT
sbatch $SCRIPT output/01_import_v1/seuratObject_IB1_BPCells_full.Rds IB1 IB1
sbatch $SCRIPT output/01_import_v1/seuratObject_RB1_BPCells_full.Rds RB1 RB1