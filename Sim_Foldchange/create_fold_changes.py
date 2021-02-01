#!/usr/bin/env python

# take on Input formatted as per Sim_Extract_IDS Script:
# transncript gene count dispersion length

# create fold_changes_r and mean_disp_r as needed for real simulation

import sys, ast, os
import random
import numpy as np	

unexpressed_genes = {}
expressed_keys = []
t_to_count = {}
t_to_disp = {}  
t_to_mean = {} 
t_to_gene = {} 

def gen_upreg():
    return np.random.choice(np.arange(2, 6, 0.5), size=1)[0]


def gen_downreg():
    return np.random.choice([ 1 / float(x) for x in np.arange(2, 6, 0.5) ], size=1)[0]

for line in open(sys.argv[1]):
     split = line.strip().split()
     t = split[0]
     g = split[1]
     count = float(split[2])
     disp = float(split[3] if split[3] != "NA" else 0)
     length = float(split[4])

     if count < 10 or disp <= 1e-4:     
         unexpressed_genes.setdefault(g,[]).append(t)
     else:
         expressed_keys.append(t)
         t_to_count[t] = count
         t_to_disp[t] = disp
         t_to_mean[t] = count * 200 / length
         t_to_gene[t] = g


fcc = int(len(expressed_keys) * 0.1) # 10 % expressed keys

fold_changed = random.sample(expressed_keys, fcc)
remaining = [x for x in expressed_keys if x not in fold_changed]

# IDs for GTF extraction

tgi = open("transcript_guide_ids_80", "w")
fold_sample = random.sample(fold_changed, int(len(fold_changed)*0.8))
for t in fold_sample:
    tgi.write("FC " + t + "\n")
remain_sample = random.sample(remaining, int(len(remaining)*0.8))
for t in remain_sample:
    tgi.write("MS " + t + "\n")

alt_transcripts_available = []
for t in fold_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

tgi1 = open("transcript_guide_ids_80_5alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
tgi2 = open("transcript_guide_ids_80_10alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")

alt_transcripts_available = []
for t in remain_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")


tgi = open("transcript_guide_ids_60", "w")
fold_sample = random.sample(fold_changed, int(len(fold_changed)*0.6))
for t in fold_sample:
    tgi.write("FC " + t + "\n")
remain_sample = random.sample(remaining, int(len(remaining)*0.6))
for t in remain_sample:
    tgi.write("MS " + t + "\n")

alt_transcripts_available = []
for t in fold_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

tgi1 = open("transcript_guide_ids_60_5alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
tgi2 = open("transcript_guide_ids_60_10alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")

alt_transcripts_available = []
for t in remain_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")


tgi = open("transcript_guide_ids_40", "w")
fold_sample = random.sample(fold_changed, int(len(fold_changed)*0.4))
for t in fold_sample:
    tgi.write("FC " + t + "\n")
remain_sample = random.sample(remaining, int(len(remaining)*0.4))
for t in remain_sample:
    tgi.write("MS " + t + "\n")

alt_transcripts_available = []
for t in fold_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

tgi1 = open("transcript_guide_ids_40_5alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
tgi2 = open("transcript_guide_ids_40_10alt", "w")
for t in random.sample(alt_transcripts_available, int(len(fold_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")

alt_transcripts_available = []
for t in remain_sample:
    g = t_to_gene[t]
    if g in unexpressed_genes:
      for at in unexpressed_genes[g]:
        alt_transcripts_available.append(at)

for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.05)):
    tgi1.write("ALT " + t + "\n")
for t in random.sample(alt_transcripts_available, int(len(remain_sample)*0.1)):
    tgi2.write("ALT " + t + "\n")


# Create Fold_Changes

fcl = open("fold_changes", "w")

t_to_cov = {}
for t in fold_changed:
    flip = random.randint(0, 1)
    if flip:
        t_to_cov[t] = gen_upreg()
    else:
        t_to_cov[t] = gen_downreg()
    fcl.write("FC " + t_to_gene[t] + " " + t + " " + str(t_to_cov[t])+ " " + str(t_to_mean[t]) + " " + str(t_to_mean[t] * t_to_cov[t]) + "\n")

 
fcr = open("fold_changes_r", "w")
mdr = open("mean_disp_r", "w")
mdr.write("Transcript Mean Disp\n")
for t in expressed_keys:

     cov = 1
     if t in t_to_cov:
         cov = t_to_cov[t]
     fcr.write( t + " " +str(cov) + "\n")

     mdr.write(t + " " + str(t_to_count[t]) + " " + str(t_to_disp[t]) + "\n")


