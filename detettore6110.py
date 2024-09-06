#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os 
import shutil
import tempfile

from lib import io
from lib import copynumbers
from lib import insertionsites
from lib import readparsing


def get_args():

    parser = argparse.ArgumentParser(
        description='A tool to infer IS6110 copy numbers and insertion sites \
        from short-read sequencing data. Default values in [].')

    parser_input = parser.add_argument_group('INPUT / OUTPUT')
    parser_settings = parser.add_argument_group('PROGRAM SETTINGS')
    
    path_to_detettore = os.path.dirname(__file__)
    
    # INPUT/OUTPUT
    parser_input.add_argument(
        'reads', nargs='+',
        help='Short reads in fasta/fastq/bam. One file for SE, two files separated by space for PE.')

    parser_input.add_argument(
        "-r", dest="ref",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0_v1.1.fasta'),
        help='Reference genome in fasta format. [resources/reference/MTBC0_v1.1.fasta]')

    parser_input.add_argument(
        '-t', dest="target",
        default=os.path.join(path_to_detettore,'resources/is_targets/IS6110.fasta'),
        help='IS consensus sequence in fasta format. [resources/is_targets/IS6110.fasta]')
    
    parser_input.add_argument(
        "-a", dest="annot",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0v1.1_PGAP_annot.gff'),
        help='Gene annotation in gff format. [resources/reference/MTBC0v1.1_PGAP_annot.gff]')
    
    # OTHER SETTINGS
    parser_input.add_argument(
        '-o', dest='outfile',
        help='Write output to this file instead of stdout.')
    
    parser_settings.add_argument(
        '-m', dest="mapq",
        type=int, default=0,
        help='Minimum mapping quality of reference-aligned reads. [0]')
    
    parser_settings.add_argument(
        '-l', dest='min_split_len',
        type=int, default=15,
        help='Minimum alignment length for splitread target hits. [10]')
    
    parser_settings.add_argument(
        '-d', dest='max_tsd_len',
        type=int, default=10,
        help='Maximum distance between left and right splitread cluster. [10]')

    parser_settings.add_argument(
        '-cl', dest="min_cl_len",
        type=int, default=10,
        help='Minimum number of anchor reads that have to cluster in order \
            to be counted in the copy number estimation [10]')

    parser_settings.add_argument(
            '-c', dest='cpus',
            type= int, default=4,
            help='Number of CPUs. [4]')

    parser_settings.add_argument(
        '--keep',
        action='store_true',
        help='Keep intermediate files.')

    args=parser.parse_args()

    return args


def main():

    args = get_args()
    
    # Mise en place ###########################################################
    
    working_dir = os.getcwd()    
    temp_dir = tempfile.mkdtemp()
    
    reads = [os.path.abspath(x) for x in args.reads]
    target = os.path.abspath(args.target)
    reference = os.path.abspath(args.ref)
    
    # Convert input bam/cram to fastq
    read_suffix = set([x.split('.')[-1] for x in reads]).pop()
    
    if len(reads) == 1 and read_suffix in ['bam', 'cram', 'sam']:

        bamfile = reads[0]        
        reads = [readparsing.bam_to_fastq(bamfile, f'{temp_dir}/reads.fastq.gz')]
        
    
    # Get reads mapping against IS6110 ########################################
    
    # Map reads against IS6110
    readparsing.mapreads(
        reads, target, 'reads_vs_IS', temp_dir, 'paf', args.cpus, k=9, m=10)

    # Extract partially mapping reads    
    partially_mapping = readparsing.get_partially_mapping(f'{temp_dir}/reads_vs_IS.paf')
    
    # Write to fastq
    anchors = readparsing.subset_fastq(partially_mapping, reads, temp_dir)


    # Estimate copy number from clustered anchor reads ########################
    anchor_clusters = copynumbers.cluster_anchors(
        [f'{temp_dir}/anchors.5.fasta', f'{temp_dir}/anchors.3.fasta'], anchors, temp_dir)
    
    copy_number = copynumbers.get_copy_number(
        anchor_clusters,  min_cluster_size = args.min_cl_len)
    
    
    # Find insertion sites in the reference genome ############################
    
    # Map partially mapping reads against reference
    readparsing.mapreads(
        [f'{temp_dir}/partially_mapping.fastq.gz'], reference, 
        'reads_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10)

    # Get split reads from reference alignment
    splitreads = readparsing.getsplitreads(
        f'{temp_dir}/reads_vs_ref.bam', temp_dir, args.min_split_len, args.mapq)


    # Map split parts against IS consensus sequences
    readparsing.mapreads(
        [f'{temp_dir}/softclipped.fasta'], target, 'softclipped_vs_IS', 
        temp_dir, 'paf', args.cpus, k=9, m=10)
    
    hits = insertionsites.parse_paf(f'{temp_dir}/softclipped_vs_IS.paf', args.min_split_len)


    # Detect clusters of split reads ##########################################

    # Cluster splitreads
    split_positions = insertionsites.cluster_splitreads(splitreads, hits)

    # Merge split clusters
    clusters = insertionsites.merge_clusters(split_positions, hits)

        
    # Write output & clean up #################################################
    io.write_output(args, clusters, copy_number, temp_dir)
    
    if args.keep:  # Copy contents of temporary to working directory
        shutil.copytree(temp_dir, os.path.join(working_dir, os.path.basename(temp_dir)))
        
    shutil.rmtree(temp_dir)


if __name__ == '__main__':
    main()

