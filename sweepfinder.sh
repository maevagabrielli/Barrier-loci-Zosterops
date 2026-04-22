#!/bin/bash

pop=$1
chr=$2
/work/project/zobo/full_genomes/sweepfinder/SF2/SweepFinder2 -lrg 50000 "$pop".SNP.frq."$chr" SFS."$pop".SNP RecRate.SNP."$chr" Out."$pop".SNP."$chr".txt
