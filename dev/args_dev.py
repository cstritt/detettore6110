#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Arguments for manual testing

Created on Tue Aug 13 09:33:42 2024

@author: cristobal
"""

class args:
    
    def __init__(self):
        
        self.reads = [
            '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5180.HS20_RL100_COV40_1.fq.gz', 
            '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5180.HS20_RL100_COV40_2.fq.gz']
        self.outfile = False
        self.ref = 'resources/reference/MTBC0_v1.1.fasta'
        self.target = 'resources/is_targets/IS6110.fasta'
        
        self.cpus = 4 
        self.mapq = 0
        self.min_split_len = 15
        self.max_tsd_len = 10
        self.min_cl_len = 10
        self.keep = False
        
args = args()

#%%
class args:
    
    def __init__(self):
        
        self.reads = [
            '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5180.HS20_RL100_COV40_1.fq.gz', 
            '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5180.HS20_RL100_COV40_2.fq.gz']
        self.outfile = False
        self.ref = 'resources/reference/H37Rv.fasta'
        self.target = 'resources/is_targets/IS6110.fasta'
        #self.annot = 'resources/reference/H37Rv.RefSeq_annot.gff'
        self.annot= False
        self.cpus = 4 
        self.mapq = 0
        self.min_split_len = 15
        self.max_tsd_len = 10
        self.min_cl_len = 10
        self.keep = False
        
args = args()

#%% bam input


#%% cigar bug
""" Easy fix: skip unmapped reads
"""

class args:
    
    def __init__(self):
        
        self.reads = ['testing/tmp_TypeError/reads.fastq.gz']
        self.outfile = False
        self.ref = 'resources/reference/MTBC0_v1.1.fasta'
        self.target = 'resources/is_targets/IS6110.fasta'
        self.annot = 'resources/reference/MTBC0_v1.1.fasta'
        self.annot= 'resources/reference/MTBC0v1.1_PGAP_annot.gff'
        self.cpus = 4 
        self.mapq = 0
        self.min_split_len = 15
        self.max_tsd_len = 10
        self.min_cl_len = 10
        self.keep = False
        
args = args()


#%% read ID bug
""" Weird problem with /1 and /2 extensions of reads

"""

class args:
    
    def __init__(self):
        
        self.reads = ['testing/tmp_KeyError/reads.fastq.gz']
        self.outfile = False
        self.ref = 'resources/reference/MTBC0_v1.1.fasta'
        self.target = 'resources/is_targets/IS6110.fasta'
        self.annot = 'resources/reference/MTBC0_v1.1.fasta'
        self.annot= 'resources/reference/MTBC0v1.1_PGAP_annot.gff'
        self.cpus = 4 
        self.mapq = 0
        self.min_split_len = 15
        self.max_tsd_len = 10
        self.min_cl_len = 10
        self.keep = False
        
args = args()
