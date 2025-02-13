#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import atexit
import os
import tempfile

from lib import io
from lib import readparsing
from lib import clusters


def get_args():

    parser = argparse.ArgumentParser(
        prog = 'dettetore6110.py',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        
        description="""
        Infer insertion sequence polymorphisms and copy numbers \
        from short-read sequencing data.\n""",
        
        epilog="""
        Example usage: """
        )

    parser_input = parser.add_argument_group('Input/Output')
    parser_settings = parser.add_argument_group('Parameters')
    path_to_detettore = os.path.dirname(__file__)
    
    # INPUT/OUTPUT
    parser_input.add_argument(
        'reads', nargs='+', 
        help='Short reads in fasta/fastq/bam. One file for SE, two files separated by space for PE.')

    parser_input.add_argument(
        '-t', dest="target", required=True,
        help='IS consensus sequence in fasta format.')
    
    parser_input.add_argument(
        '-o', dest='outpath',
        help='Path to output directory.')
    
    parser_input.add_argument(
        '-p', dest='prefix',
        help='Prefix for output files.')
    
    parser_input.add_argument(
        "-r", dest="reference",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0_v1.1.fasta'),
        help='Reference genome in fasta format.')

    parser_input.add_argument(
        "-a", dest="annot",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0v1.1_PGAP_annot.gff'),
        help='Gene annotation in gff format.')
    
    
    # OTHER SETTINGS
    parser_settings.add_argument(
        '-al', dest='min_anchor_len',
        type=int, default=20,
        help='Minimum length of the read part that maps outside the IS.')

    parser_settings.add_argument(
        '-hl', dest='min_hit_len',
        type=int, default=20,
        help='Minimum length of the read part that maps to the IS.')
    
    parser_settings.add_argument(
        '-cs', dest='min_cluster_size',
        type=int, default=5,
        help='Minimum number of anchor reads in a cluster.')
    
    parser_settings.add_argument(
        '-tsd', dest='tsd_len',
        nargs='+', type=int, default=[3,4],
        help='Alowable length of the target site duplication.')
    
    parser_settings.add_argument(
        '-c', dest='cpus',
        type= int, default=4,
        help='Number of CPUs.')
    
    parser_settings.add_argument(
        '--keep', dest='pref', type=str, default=False,
        help='Keep intermediate files in folder <pref>_tmp.')

    args=parser.parse_args()

    return args


def main():
    
    args = get_args()
    
    # Le mise-en-place ########################################################
    working_dir = os.getcwd()    
    reads = [os.path.abspath(x) for x in args.reads]
    target = os.path.abspath(args.target)
    temp_dir = tempfile.mkdtemp()
    atexit.register(io.exit_handler, args, temp_dir, working_dir)

    # Convert input bam/cram to fastq
    read_suffix = set([x.split('.')[-1] for x in reads]).pop()

    if len(reads) == 1 and read_suffix in ['bam', 'cram', 'sam']:
        bamfile = reads[0]        
        reads = [readparsing.bam_to_fastq(bamfile, f'{temp_dir}/reads.fastq.gz')]


    # Map reads against IS target ###########################################
    readparsing.mapreads(
        reads, target, 'reads_vs_IS', temp_dir, 'paf', args.cpus, k=9, m=10
    )

    # Create read dictionary
    read_d = readparsing.parse_paf(f'{temp_dir}/reads_vs_IS.paf')

    # Add anchor and hit parts of the reads to read dictionary
    readparsing.add_seqs_to_read_dict(read_d, reads, temp_dir)

    # Cluster anchors based on exact identity of anchor part adjoining IS
    anchor_clusters = clusters.cluster_anchors(read_d, args.min_anchor_len)

    # Create cluster dictionary and anchor consensus sequences
    cluster_d = clusters.parse_clusters(read_d, anchor_clusters, temp_dir, args)


    # Identify reference positions ##########################################
    readparsing.mapreads(
        [f'{temp_dir}/anchor_consensi.fasta'], args.reference, 'reads_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10
        )

    overlaps = clusters.find_overlaps(f'{temp_dir}/reads_vs_ref.bam', args.tsd_len)
    clusters.add_ref_coordinates_to_clusters(cluster_d, f'{temp_dir}/reads_vs_ref.bam')


    # Write output ##########################################################
    io.write_cluster_output(cluster_d, args)
    io.write_reference_output(overlaps, cluster_d, args)
    
    
if __name__ == '__main__':
    main()
