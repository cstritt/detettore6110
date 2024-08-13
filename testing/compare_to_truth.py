#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Compare detettore output to truth (liftover of IS annotation with 
odgi position)


Created on Wed Apr 24 14:03:55 2024

@author: cristobal
"""


import pandas

truth_path = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/mtbc30.is6110_sites.Anc01_liftover.tsv'

truth = pandas.read_csv(truth_path, sep='\t')
