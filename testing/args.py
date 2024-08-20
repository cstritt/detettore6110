#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Arguments for manual testing

Created on Tue Aug 13 09:33:42 2024

@author: cristobal
"""

class args:
    
    def __init__(self):
        
        self.fastq = ['/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5079.HS25_RL150_COV50_1.fq.gz', 
                      '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5079.HS25_RL150_COV50_2.fq.gz']
        self.outfile = False
        self.ref = 'resources/reference/MTBC0_v1.1.fasta'
        self.target = 'resources/is_targets/IS6110.fasta'
        
        self.cpus = 4 
        self.mapq = 0
        self.min_split_len = 15
        self.max_tsd_len = 10
        self.min_cl_len = 10
        self.keep = True
        
args = args()