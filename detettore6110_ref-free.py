#!/usr/bin/env python3

# %%

import argparse
import os
import shutil
import tempfile

from lib import io
from lib import readparsing
from lib import clusters


def get_args():

    parser = argparse.ArgumentParser(
        description='Infer insertion sequence polymorphisms from \
            short-read sequencing data. Default values in [].')

    parser_input = parser.add_argument_group('INPUT / OUTPUT')
    parser_settings = parser.add_argument_group('PROGRAM SETTINGS')
    

    path_to_detettore = os.path.dirname(__file__)
    
    # INPUT/OUTPUT
    parser_input.add_argument(
        'reads', nargs='+', required=True,
        help='Short reads in fasta/fastq/bam. One file for SE, two files separated by space for PE.')

    parser_input.add_argument(
        '-t', dest="target", required=True,
        help='IS consensus sequence in fasta format.')


    # Optional reference
    parser_input.add_argument(
        "-r", dest="ref",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0_v1.1.fasta'),
        help='Reference genome in fasta format. [resources/reference/MTBC0_v1.1.fasta]')

    parser_input.add_argument(
        "-a", dest="annot",
        default=os.path.join(path_to_detettore, 'resources/reference/MTBC0v1.1_PGAP_annot.gff'),
        help='Gene annotation in gff format. [resources/reference/MTBC0v1.1_PGAP_annot.gff]')
    
    
    # OTHER SETTINGS
    parser_settings.add_argument(
        '-al', dest='min_anchor_len',
        type=int, default=20,
        help='Minimum length of the read part that maps outside the IS. [20]')

    parser_settings.add_argument(
        '-hl', dest='min_hit_len',
        type=int, default=20,
        help='Minimum length of the read part that maps to the IS. [20]')
    
    parser_settings.add_argument(
        '-cs', dest='min_cluster_size',
        type=int, default=5,
        help='Minimum number of anchor reads in a cluster. [5]')

    parser_input.add_argument(
        '-o', dest='outfile',
        help='Write output to this file instead of stdout.')

    parser_settings.add_argument(
        '--keep', dest='pref', type=str, default=False,
        help='Keep intermediate files in folder <pref>_tmp.')
    
    parser_settings.add_argument(
        '-c', dest='cpus',
        type= int, default=4,
        help='Number of CPUs. [4]')

    args=parser.parse_args()

    return args


#args = get_args()
#%%
class args:
    def __init__(self):
        self.reads = ['testing/some_reads.fastq.gz']
        self.target = 'resources/is_targets/IS6110.fasta'
        self.cpus=4
        self.min_anchor_len=20
        self.min_hit_len=20
        self.min_cluster_size=5
        self.tsd_len = [3,4]
        
        self.reference = 'resources/reference/MTBC0_v1.1.fasta'
        self.annot = 'resources/reference/MTBC0v1.1_PGAP_annot.gff'
        
args = args()

# Le mise-en-place ########################################################

working_dir = os.getcwd()    
reads = [os.path.abspath(x) for x in args.reads]
target = os.path.abspath(args.target)
#temp_dir = tempfile.mkdtemp()
#atexit.register(io.exit_handler, args, temp_dir, working_dir)
temp_dir = 'testing/tmp'
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
os.mkdir(temp_dir)

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

# Add anchor and hit parts of the reads
readparsing.add_seqs_to_read_dict(read_d, reads, temp_dir)

# Cluster anchors based on exact identity of first 20/50 bp
anchor_clusters = clusters.cluster_anchors(read_d, 50)

# Create cluster dictionary and consensus sequences
cluster_d = clusters.parse_clusters(read_d, anchor_clusters, temp_dir, args)

out = clusters.write_results(cluster_d, args, temp_dir)

# Write cluster output




# Identify reference positions ##########################################
readparsing.mapreads(
    [f'{temp_dir}/anchor_consensi.fasta'], args.reference, 'reads_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10
    )

overlaps = clusters.find_overlaps(f'{temp_dir}/reads_vs_ref.bam', args.tsd_len)

io.write_reference_output(overlaps, args, temp_dir)
