#!/bin/bash
#SBATCH -e projection_all.sh.e
#SBATCH -o projection_all.sh.o

# input needs to be tab separated
sed -e 's/,/\t/g' scans_vir_pall_pi_dxy_fst_all.csv >scans_vir_pall_pi_dxy_fst_all.tab

# project onto the zebra finch pseudo-chromosomes with a perl script available at: https://osf.io/uw6mb/ 
perl /work/project/zobo/full_genomes/scans_diversity/ProjectionPositionsPseudoMolecule_windows.pl -a /work/project/zobo/full_genomes/scans_diversity/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_Zosterops_borbonicus_ZeFi-ZoBo+DeCoSTAR.agp -i scans_vir_pall_pi_dxy_fst_all.tab

grep "chromosome" scans_vir_pall_pi_dxy_fst_all.tab.pseudoK >scans_vir_pall_pi_dxy_fst_all.tab.chr.pseudoK

head -n 1 scans_vir_pall_pi_dxy_fst_all.tab >header
echo chr >chr
echo pos >pos
paste header chr pos >header2
cat header2 scans_vir_pall_pi_dxy_fst_all.tab.chr.pseudoK >scans_vir_pall_pi_dxy_fst_all.tab.chr.header.pseudoK

rm header header2
