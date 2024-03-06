#!/bin/bash
#SBATCH --job-name=testing_ei_pipeline
#SBATCH --output=sm_test.out
#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=isaac.wade@icr.ac.uk
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

. "/opt/software/applications/anaconda/3/etc/profile.d/conda.sh" && conda activate /data/scratch/DGE/DUDGE/PREDIGEN/iwade/SMARCA4_related/workflow/envs/SMARCA4_related

CONFIG_FILE="/data/scratch/DGE/DUDGE/PREDIGEN/iwade/EI_testing/config/config.yaml"

# Prepare Dir Structure
mkdir -p test_slurm

[ -e pipeline_complete.txt ] && rm pipeline_complete.txt

snakemake --snakefile ../../workflow/rules/find_var_carriers.smk \
	--configfile "$CONFIG_FILE" \
	--cores 1 \
	--jobs 50 \
	--default-resources "time='24:00:00'" \
	--latency-wait 36000 \
	-n

# Tidy .out files
if ls slurm*.out 1>/dev/null 2>&1; then mv slurm*.out test_slurm/; fi

