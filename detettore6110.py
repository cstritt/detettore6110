#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" detettore6110 

A tool to detect IS6110 polymorphisms in the Mycobacterium tuberculosis complex 
based on split reads 

To do:
    - copy number output
    - merging of L and R clusters: large negative values
    - 


"""


#%% Libraries

import argparse
import os 
import shutil
import sys

from scripts import insertion_sites_module
from scripts import copy_numbers_module


#%% Parse arguments

def get_args():

    parser = argparse.ArgumentParser(
        description='Detect IS6110 polymorphisms from short-read data. \
            Default values in square brackets.')


    parser_input = parser.add_argument_group('INPUT / OUTPUT')
    parser_settings = parser.add_argument_group('PROGRAM SETTINGS')
    
    path_to_detettore = os.path.dirname(__file__)

    # INPUT/OUTPUT
    parser_input.add_argument(
        '-f', dest='FASTQ',
        required=True,
        nargs = '+',
        help='Short reads in fastq format. One file for SE, two files separated by space for PE.')
    
    parser_input.add_argument(
        '-o', dest='OUTPREF',
        required=True,
        help='Prefix for output files.')

    parser_input.add_argument(
        "-r", dest="REF",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0_v1.1.fasta'),
        help='Reference genome in fasta format. [resources/reference/MTBC0_v1.1.fasta]')

    parser_input.add_argument(
        '-t', dest="TARGETS",
        default=os.path.join(path_to_detettore,'resources/is_targets/IS6110.fasta'),
        help='IS consensus sequences in fasta format. [resources/is_targets/IS6110.fasta]')

    parser_input.add_argument(
        '-a', dest="ANNOT",
        default=os.path.join(path_to_detettore,'resources/is_annotation/MTBC0.ISEScan.gff'),
        help='IS annotation in gff format. [resources/is_annotation/MTBC0.ISEScan.gff]')
    
    parser_input.add_argument(
        '-g', dest="GENES",
        help='Reference gene annotation in gff format. [resources/reference/MTBC0v1.1_PGAP_annot.gff]',
        default=os.path.join(path_to_detettore,'resources/reference/MTBC0v1.1_PGAP_annot.gff'))
    
    parser_input.add_argument(
        '--to-file', dest='to_file',
        action=('store_true'),
        help='Write to file rather than stdout')
    
    # OTHER SETTINGS
    parser_settings.add_argument(
            '-c', dest='CPUS',
            type= int, default=4,
            help='Number of CPUs. [4]')
    
    parser_settings.add_argument(
        '-m', dest="MAPQ",
        type=int, default=0,
        help='Minimum mapping quality of reference-aligned reads. [0]')
    
    parser_settings.add_argument(
        '-l', dest='MIN_LEN',
        type=int, default=15,
        help='Minimum alignment length for splitread target hits. [15]')
    
    parser_settings.add_argument(
        '-d', dest='MAX_LR_DIST',
        type=int, default=10000,
        help='Maximum distance between left and right splitread cluster. [10000]')

    parser_settings.add_argument(
        '--keep',
        action='store_true',
        help='Keep intermediate files.')

    args=parser.parse_args()

    return args


#%% MAIN

def main():

    args = get_args()
    params = insertion_sites_module.mise_en_place(args)
    
    if os.path.exists(params.tmp):
        shutil.rmtree(params.tmp)
    os.mkdir(params.tmp)
    

    #%% Get reads mapping against IS6110
    
    # Map reads against IS6110
    #sys.stderr.write('Mapping reads against IS ...\n')
    insertion_sites_module.mapreads(params.FASTQ, params.TARGETS, 'reads_vs_IS', 4, 'paf', params.tmp, k=9, m=5)

    # Extract partially mapping reads    
    partially_mapping = insertion_sites_module.get_partially_mapping('{}/reads_vs_IS.paf'.format(params.tmp), params.FASTQ, args)
    
    # Write to fastq
    anchors = insertion_sites_module.subset_fastq(partially_mapping, params.FASTQ, params)


    #%% Estimate copy number from clustered anchor reads
    #sys.stderr.write('Estimating copy number from clustered anchor reads ...\n')
    
    clusters = copy_numbers_module.cluster_anchors(
        ['{}/anchors.5.fasta'.format(params.tmp), '{}/anchors.3.fasta'.format(params.tmp)], anchors, params)
    
    copy_number = copy_numbers_module.get_copy_number(clusters, params, min_cluster_size = 10)
    
    with open(params.OUTPREF + '.copy_number.txt', 'w') as f:
        f.write(str(copy_number))


    #%% Find insertion sites in the reference genome
    #sys.stderr.write('Detecting insertion sites in the reference genome ...\n')
    
    # Map partial hits against reference
    insertion_sites_module.mapreads(['{}/partially_mapping.fastq.gz'.format(params.tmp)], params.REF, '{}/{}'.format(params.tmp,params.OUTPREF), params.CPUS, 'sam', params.tmp)

    # Get split reads from reference alignment
    splitreads = insertion_sites_module.getsplitreads('{}/{}.bam'.format(params.tmp,params.OUTPREF), params.tmp, params.MIN_LEN, params.MAPQ)

    # Map split parts against IS consensus sequences
    hits = insertion_sites_module.ISmap('{}/softclipped.fasta'.format(params.tmp), params.TARGETS, '{}/{}'.format(params.tmp,params.OUTPREF), min_aln_len = 10, k = 9, w = 5)


    #%% Detect clusters of split reads

    # Cluster splitreads
    split_positions = insertion_sites_module.cluster_splitreads(splitreads, hits)

    # Merge split clusters
    clusters = insertion_sites_module.merge_clusters(split_positions, hits)

        
    #%% Write output & clean up
    insertion_sites_module.write_output(params, clusters, params.tmp, to_file=args.to_file) 
    
    if not args.keep:
        shutil.rmtree(params.tmp)


#%% 
if __name__ == '__main__':
    main()


