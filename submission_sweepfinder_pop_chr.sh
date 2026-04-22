#!/bin/bash

for pop in HN HS NE S W
do
cat chr.list | while read chr
do
sbatch -J SF."$pop"."$chr" -e SF."$pop"."$chr".e -o SF."$pop"."$chr".o sweepfinder.sh $pop $chr
done
done
