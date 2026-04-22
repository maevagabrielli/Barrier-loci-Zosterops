#!/bin/bash
#SBATCH -o sliding_window_SNP.sh.o
#SBATCH -e sliding_window_SNP.sh.e
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16

module load devel/python/Python-2.7.18
# input
python /work/project/zobo/full_genomes/scans_diversity/parseVCF.py -i /home/mgabrielli/work/reviewZosteropsEvolution/scans_diversity_noF/vir_pall_merged_joint_bwa_mem_mdup.filtered.noindel.PASS.DP5-150.nomiss.SNP.noF.vcf.gz -o vir.pall.filtered.noindel.PASS.DP5-150.nomiss.SNP.noF.geno
#lasted 30minutes
gzip vir.pall.filtered.noindel.PASS.DP5-150.nomiss.SNP.noF.geno

# use popgenWindows by Simon Martin, in 50kb windows
python /work/project/zobo/full_genomes/scans_diversity/popgenWindows.py -g vir.pall.filtered.noindel.PASS.DP5-150.nomiss.SNP.noF.geno.gz --windType coordinate -w 50000 -o scans_vir_pall_pi_dxy_fst_all.csv -f phased --popsFile pops.vir.pall.txt -p VIR -p PALL --ploidy 2 -T 16
