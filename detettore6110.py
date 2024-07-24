#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" detettore6110 

A tool to detect IS6110 polymorphisms in the Mycobacterium tuberculosis complex 
based on split reads 


"""


#%% Libraries

import argparse
import os 
import shutil

from scripts import strumenti


#%% Parse arguments

def get_args():

    parser = argparse.ArgumentParser(
        description='Detect IS6110 polymorphisms from short-read data')


    parser_input = parser.add_argument_group('INPUT / OUTPUT')
    parser_settings = parser.add_argument_group('PROGRAM SETTINGS')

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
        default='resources/reference/MTBC0_v1.1.fasta',
        help='Reference genome in fasta format.')

    parser_input.add_argument(
        '-t', dest="TARGETS",
        default='resources/is_targets/IS6110.fasta',
        help='IS consensus sequences in fasta format.')

    parser_input.add_argument(
        '-a', dest="ANNOT",
        default='resources/is_annotation/Anc01.fasta.MTBC0_renamed.gff',
        help='IS annotation in gff format.')
    
    parser_input.add_argument(
        '-g', dest="GENES",
        help='Reference gene annotation in gff format.',
        default='resources/reference/MTBC0v1.1_PGAP_annot.gff')
    
    parser_input.add_argument(
        '--to-file',
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
        '--keep',
        action='store_true',
        help='Keep intermediate files.')

    args=parser.parse_args()

    return args


#%% MAIN

def main():

    args = get_args()
    params = strumenti.mise_en_place(args)
    
    if os.path.exists(params.tmp):
        shutil.rmtree(params.tmp)
    os.mkdir(params.tmp)
    

    #%% Get reads mapping against IS6110
    
    # Map reads against IS6110
    strumenti.mapreads(params.FASTQ, params.TARGETS, 'reads_vs_IS', 4, 'paf', params.tmp, k=9, m=5)

    # Extract partially mapping reads    
    partially_mapping = strumenti.get_partially_mapping(f'{params.tmp}/reads_vs_IS.paf', params.FASTQ, args)
    
    # Write to fastq
    anchors = strumenti.subset_fastq(partially_mapping, params.FASTQ, params)


    #%% Estimate copy number from clustered anchor reads
    clusters = strumenti.cluster_anchors(
        [f'{params.tmp}/anchors.5.fasta', f'{params.tmp}/anchors.3.fasta'], anchors, params)
    
    
    #%% Find insertion sites in the reference genome

    # Map partial hits against reference
    strumenti.mapreads([f'{params.tmp}/partially_mapping.fastq.gz'], params.REF, f'{params.tmp}/{params.OUTPREF}', params.CPUS, 'sam', params.tmp)

    # Get split reads from reference alignment
    splitreads = strumenti.getsplitreads(f'{params.tmp}/{params.OUTPREF}.bam', params.tmp, params.MIN_LEN, params.MAPQ)

    # Map split parts against IS consensus sequences
    hits = strumenti.ISmap(f'{params.tmp}/softclipped.fasta', params.TARGETS, f'{params.tmp}/{params.OUTPREF}', min_aln_len = 10, k = 9, w = 5)


    #%% Detect clusters of split reads

    # Cluster splitreads
    split_positions = strumenti.cluster_splitreads(splitreads, hits)

    # Merge split clusters
    clusters = strumenti.merge_clusters(split_positions, hits)

        
    #%% Write output & clean up
    strumenti.write_output(params, clusters, params.tmp, to_file=False) 
    
    if not args.keep:
        shutil.rmtree(params.tmp)


#%% 
if __name__ == '__main__':
    main()


