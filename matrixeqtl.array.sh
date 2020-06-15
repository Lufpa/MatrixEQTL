#!/bin/bash
#SBATCH --mem=20000
#SBATCH --time=24:00:00 --qos=1day
#SBATCH --job-name=eqtsub
#SBATCH --cpus-per-task=1   #make sure to modify $cpus too!!!
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --array=4

##Requires:
# a sampleinfo sheet with 7 tab separated columns 
# covariatfile, analysis(g or gxe), vstfile, normalized (yes or no), #PCs from kinship, output filehandle, run without covariates (when using a batch adjusted count matrix)

# it should be run from the folder than contains the sampleinfo file as well as all others.

# change manually the genotype file from matrixeqtl.subsets.R
# make sure the VCF is in BED format (use plint --vcf --make-bed --out --allow-extra-chr)
# and change bim file to recode chromosomes as numeric for gtca and/or ldka

# Super important: don't forget to conda activate R_MatrixEQTL environment before launching the job

set -e
date
echo "starting matrixeqtl for subset of samples"

info=infosheet.txt

source ~/scripts/MatrixEqtl/matrixeqtl.launchR.sh $info

echo "done!"
date

