#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import atexit
import os
import shutil
import tempfile

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
        Example usage:
        
        detettore6110.py testing/some_reads.fastq.gz \
            -t resources/is_targets/IS6110.fasta \
            -o testing/results \
            -p some_reads
            
        With reference genome and annotation:
        
        detettore6110.py testing/some_reads.fastq.gz \
            -t resources/is_targets/IS6110.fasta \
            -r resources/reference/MTBC0_v1.1.fasta \
            -a resources/reference/MTBC0v1.1_PGAP_annot.gff \
            -o testing/results \
            -p some_reads
        
        """
        )

    parser_input = parser.add_argument_group('Input/Output')
    parser_settings = parser.add_argument_group('Parameters')
    path_to_detettore = os.path.dirname(__file__)
    
    # Input/Output
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
        help='Reference genome in fasta format.')

    parser_input.add_argument(
        "-a", dest="annot",
        help='Gene annotation in gff format.')
    
    # Parameters
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
        type=int, default=3,
        help='Minimum number of anchor reads in a cluster.')
    
    parser_settings.add_argument(
        '-k', dest='seed_len',
        type=int, default=20,
        help='Require k exact matches next to the breakpoint for anchor reads to cluster.')
    
    parser_settings.add_argument(
        '-tsd', dest='tsd_len',
        nargs='+', type=int, default=[3,4],
        help='Alowable length of the target site duplication.')
    
    # Other settings
    parser_settings.add_argument(
        '-c', dest='cpus',
        type= int, default=4,
        help='Number of CPUs.')
    
    parser_settings.add_argument(
        '--keep', action = 'store_true',
        help='Keep intermediate files in folder <pref>_tmp.')

    args=parser.parse_args()

    return args


def exit_handler(args, temp_dir):
    """ Cleanup after program finish. If --keep is given, copy contents of 
    temporary directory to working directory before deleting it"""
    
    if args.keep:  # Copy contents of temporary to output directory
        dest = os.path.join(args.outpath, args.prefix + '_intermediate_files')
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(temp_dir, dest )
    shutil.rmtree(temp_dir)


def main():
    
    args = get_args()
    
    # Le mise-en-place ######################################################
    temp_dir = tempfile.mkdtemp()
    atexit.register(exit_handler, args, temp_dir)
    target = os.path.abspath(args.target)
    reads = readparsing.Reads(args, temp_dir)
    
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    
    # Map reads against IS target ###########################################
    readparsing.mapreads(reads.fastq, target, 'reads_vs_IS', temp_dir, 'paf', args.cpus, k=9, m=10)  # Map reads against target IS
    reads.parse_paf(f'{temp_dir}/reads_vs_IS.paf', args.min_anchor_len, args.min_hit_len)  # Extract reads that reach into the IS
    reads.add_seqs_to_read_dict(reads.fastq, temp_dir)  # Re-traverse reads and extract the anchor sequences

    # Identify anchor clusters ##############################################
    anchor_clusters = clusters.Clusters(args, temp_dir)
    anchor_clusters.cluster_anchors(reads, args.seed_len)  # Cluster anchor sequences based on exact identity of anchor part adjoining IS
    anchor_clusters.parse_clusters(reads)  # Align reads and get anchor consensus sequences

    # Identify reference positions ##########################################
    if args.reference:
        readparsing.mapreads([f'{temp_dir}/anchor_consensi.fasta'], args.reference, 'anchors_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10)  # Map anchors against reference
        readparsing.mapreads(reads.fastq, args.reference, 'reads_vs_ref', temp_dir, 'bam', args.cpus)  # Map all reads against reference

        anchor_clusters.add_ref_coordinates_to_clusters(f'{temp_dir}/anchors_vs_ref.bam')
        overlaps = clusters.find_overlaps(f'{temp_dir}/anchors_vs_ref.bam', args.tsd_len)
    
    # Write output ##########################################################
    anchor_clusters.write_cluster_output(args)
    if args.reference:
        anchor_clusters.write_reference_output(overlaps, f'{temp_dir}/reads_vs_ref.bam', args)
    
if __name__ == '__main__':
    main()
