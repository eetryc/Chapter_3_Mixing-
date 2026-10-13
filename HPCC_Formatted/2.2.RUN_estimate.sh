#!/bin/bash
#SBATCH --job-name=ckmr_estimation
#SBATCH --output=ckmr_estimation_%j.out
#SBATCH --error=ckmr_estimation_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trycemil@msu.edu

# Load any required R or compilation modules
module purge
module load R/4.3.2-foss-2023a 2>/dev/null || module load R/4.3.2
module load JAGS/4.3.2-foss-2023a 2>/dev/null || module load JAGS 2>/dev/null || true

# Explicitly set your user library path
export R_LIBS_USER="/mnt/home/trycemil/R/x86_64-pc-linux-gnu-library/4.3"

echo "=== Starting CKMR Simulation & Estimation ==="
date

# Run the updated R script
Rscript run_estimation.R

echo "=== Estimation Finished Successfully ==="
date
