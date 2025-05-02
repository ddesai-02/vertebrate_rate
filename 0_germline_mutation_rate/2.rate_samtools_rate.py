# -*- coding: utf-8 -*-
"""
This script ccalculate a mutation rate.

Lucie Bergeron
    Devan Desai
03.12.18
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Import species names:
species = pd.read_csv('sp_here.txt', sep=' ', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    # Directories:
    direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
    direct_handling="{}/vcf_handling/".format(path, sp)
    direct="{}/{}/".format(path, sp)
    # Dictionary:
    f = open('/{}/pedigree.ped'.format(direct))
    trio_dir = {}
    for line in f:
        off = line.split()[1]
        fa = line.split()[2]
        mo = line.split()[3]
        name = off
        if name not in trio_dir:
            trio_dir[name] = []
        trio_dir[name].append((off, fa, mo))

    # Import nb mutation and fnr
    nb_denovo = pd.read_csv('{}/denovo.txt'.format(sp_target),sep=' ', index_col=None, header=None)
    fnr = pd.read_csv('{}fnr.txt'.format(direct_denovo),sep=' ', index_col=None, header=None)
    call = pd.read_csv('{}callability.txt'.format(direct_denovo),sep=' ', index_col=None, header=None)
    # Find the overall fnr:
    fnr_all=round(sum(fnr[2])/sum(fnr[1]),5)
    # Alpha allelic balance AND site filter
    a_RP=0.0025696749
    a_MQRS=0.0331547715
    a_FS=0.0182736502
    alpha_all = 1 - ((1-a_FS)*(1-a_RP)*(1-a_MQRS)*(1-fnr_all))
    print('alpha={}'.format(alpha_all))
    # Find the mutation rate per trios:
    file = open('{}/mutation_rate.txt'.format(sp_target),'w')
    for name in trio_dir:
        nb_mut=nb_denovo.loc[nb_denovo[0] ==name][1]
        C=call.loc[call[0] ==name][1]
        alpha=alpha_all
        mut_rate = (float(nb_mut))/(2*float(C)*(1-float(alpha)))
        file.write('{} {} \n'.format(name, mut_rate))
    file.close()
