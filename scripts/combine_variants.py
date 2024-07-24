#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Combine variants obtained with detettore TB into a single table, adding 
gene information 


Created on Fri Jun 28 09:48:00 2024
@author: cristobal


Input:
    - tab-sep list of detettore result files
    - reference gene annotation in gff
    
    - apply filters at this stage?
    
    
Output:
    - vcf-like table with all samples and gene context


        #'left_gene',
        #'left_strand',
        #'left_distance',
        #'left_description',
     
        #'right_gene',
        #'right_strand',
        #'right_distance',
        #'right_description'


"""

