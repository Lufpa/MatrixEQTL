#/bin/bash

infofile=$1

covfile=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $1 }' $infofile`
analysis=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $2 }' $infofile`
vst=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $3 }' $infofile`
nor=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $4 }' $infofile`
PCs=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $5 }' $infofile`
file=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $6 }' $infofile`
covfree=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $7 }' $infofile`

#Rscript ~/scripts/MatrixEqtl/matrixeqtl.subsets.R $covfile $analysis $vst $nor $PCs $file $covfree

	

	#if ruuning matrixeQTL for multi-tissue analylis the source code is different

	Rscript ~/scripts/MatrixEqtl/matrixeqtl.multitissue.R $covfile $analysis $vst $nor $PCs $file $covfree
