#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Arguments for manual testing

Created on Tue Aug 13 09:33:42 2024

@author: cristobal
"""

class args:
    
    def __init__(self):
        
        self.FASTQ = ['/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5079.HS25_RL150_COV50_1.fq.gz', 
                      '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/reads/CCDC5079.HS25_RL150_COV50_2.fq.gz']
        self.OUTPREF = 'CCDC5079.HS25_RL150_COV50'
        self.REF = 'resources/reference/MTBC0_v1.1.fasta'
        self.TARGETS = 'resources/is_targets/IS6110.fasta'
        self.ANNOT = 'resources/is_annotation/MTBC0.ISEScan.gff'
        self.GENES = 'resources/reference/MTBC0v1.1_PGAP_annot.gff'
        self.to_file = False
        
        self.CPUS = 4 
        self.MAPQ = 0
        self.MIN_LEN = 15
        self.MAX_LR_DIST = 10000
        self.keep = True
        
args = args()