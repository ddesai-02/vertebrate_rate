# -*- coding: utf-8 -*-
"""
This script calls variants with mpilup for each
Lucie Bergeron
    Devan Desai
03.12.18
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Import species names:
species = pd.read_csv('sp_here.txt', sep=' ', index_col=None, header=None)
pedigree = pd.read_csv('pedigree.txt', sep='\t', index_col=None, header=None)
genome = pd.read_csv('ref.txt', sep='\t', index_col=None, header=None)

################
# What you run #
################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    genome_target = genome.loc[genome[0] == sp_target][1].to_list()[0]
    # Create directory
    mkdir = "mkdir {}".format(sp)
    subprocess.call(mkdir, shell=True)
    # Directories:
    direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
    direct_bam="{}/{}/bam_files/".format(path, sp)
    direct_ref="{}/{}/ref_fasta/".format(path, sp)
    direct = "{}/{}/".format(path, sp)
    # Dictionary:
    pedigree_sp=pd.read_csv('{}/pedigree.ped'.format(direct), sep=' ', index_col=None, header=None)
    for sample in range(0,len(pedigree_sp)):
        fa = pedigree_sp.iloc[sample,2]
        mo = pedigree_sp.iloc[sample,3]
        off = pedigree_sp.iloc[sample,1]
        denovo_to_check=pd.read_csv('{}/data_denovo_{}.tab'.format(direct_denovo, off), sep='\t', index_col=None)
        file = open('{}/{}_samtools.sh'.format(sp_target, off),'w')
        file.write('#!/bin/bash \n')
        file.write('#SBATCH --account={} \n'.format(account))
        file.write('#SBATCH --mem 10G \n')
        file.write('#SBATCH --cpus-per-task=1 \n')
        file.write('#SBATCH --time=11:59:00 \n')
        for line in range(0,len(denovo_to_check)):
            chrom=denovo_to_check.iloc[line,0]
            pos=denovo_to_check.iloc[line,1]
            file.write('samtools mpileup -ugf {}{} -r {}:{}-{} {}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam | bcftools call -m | tail -1 | cut -f1,2,4,5,10,11,12 >> {}/{}/de_novo_mutation/{}_samtools.txt \n'.format(direct_ref, genome_target, chrom, pos, pos, direct_bam, fa, direct_bam, mo, direct_bam, off, direct, sp, off))
        file.close()
        sub_cmd = "sbatch -o {}/{}_samtools.out {}/{}_samtools.sh".format(sp, off, sp, off)
        subprocess.call(sub_cmd, shell=True)
