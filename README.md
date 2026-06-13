# Barrier-loci-Zosterops
Scripts for the manuscript under review "Accounting for recombination rate variation improves inference of barrier loci and reveals the role of both natural and sexual selection in an incipient bird radiation"

## 0) recombination rate 
src_rec_rate_genome.R

## 1) Fst sans
sliding_window_SNP.sh
projection.sh

## 2) Outlier Fst accounting for recombination rate, and Manhattan plot
src_manhattan_plots_fst_recombination.R

## 3) sweepfinder
submission_sweepfinder_pop_chr.sh
sweepfinder.sh

## 4) DILS
Run_DILS.sh
config.yaml (example file)

## 5) Statistical test overlap outliers
circular_bootstrap_trio.R
circular_bootstrap_pairwise.R
