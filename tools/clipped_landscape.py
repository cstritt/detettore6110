#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" 

Created on Wed Apr 12 13:46:49 2023
@author: cristobal
"""

import pysam

from scripts import strumenti


bamfile = "/home/cristobal/TB/analyses/detettore/SE_100/G66401.bam"


#%%


split_positions = {
    "left" : {},
    "right" : {}
    }

pybam = pysam.AlignmentFile(bamfile, "rb")

for read in pybam.fetch():
    
    cigar = read.cigartuples
    
    if len(cigar) > 1: 
    
        first, last = cigar[0], cigar[-1]
        
        if first[0] in (4,5):
            pos = read.reference_start
            side = "right"
            
        elif last[0] in(4,5):
            pos = read.reference_end
            side = "left"
            
        if pos not in split_positions[side]:
            split_positions[side][pos] = 0
                
        split_positions[side][pos] += 1
            
pybam.close()


remove = {}

for side in ("left", "right"):
    
    remove[side] = []
    
    for pos in split_positions[side]:
        if split_positions[side][pos] == 1:
            remove[side].append(pos)
            
    for pos in remove[side]:
        split_positions[side].pop(pos, None)
        
        
        
        
 
#%%


with open("/home/cristobal/TB/projects/IS_GWAS/analyses/detettore/cliplandscape/G66401.cliplandscape.tsv", "w") as f:
    for side in split_positions:
        for pos in split_positions[side]:
            outline = [str(pos), str(side), str(split_positions[side][pos])]
            f.write('\t'.join(outline) + '\n')







#%% WHat happens with the hardclipped part of reads?

import pysam
bamfile = "/home/cristobal/TB/projects/IS_GWAS/analyses/detettore/G46697.bam"

hardclipped = {}


pybam = pysam.AlignmentFile(bamfile, "rb")

for read in pybam.fetch():
    
    cigar = read.cigartuples
    
    if len(cigar) > 1: 
    
        first, last = cigar[0], cigar[-1]
        
        if first[0] == 5 or last[0] == 5:
            
            if read.qname not in hardclipped:
                hardclipped[read.qname] = []
                
            hardclipped[read.qname].append(read)
            
    
pybam.close()

#%%


