#!/usr/bin/env python3
#%%
import argparse
import os
import tempfile

from lib import classi_e_funzioni as cef
from lib import readparsing
from lib import reads
from lib import clusters

#%% MISE EN PLACE
class args:
    def __init__(self):
        
        self.reads = ['/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/v2/G22/51/3/G22513.cram']
        self.target = 'resources/is_targets/IS6110.fasta'
        self.cpus = '4'
        

args = args()

working_dir = os.getcwd()    
temp_dir = tempfile.mkdtemp()

reads = [os.path.abspath(x) for x in args.reads]
target = os.path.abspath(args.target)


# Convert input bam/cram to fastq
read_suffix = set([x.split('.')[-1] for x in reads]).pop()

if len(reads) == 1 and read_suffix in ['bam', 'cram', 'sam']:
    bamfile = reads[0]        
    reads = [readparsing.bam_to_fastq(bamfile, f'{temp_dir}/reads.fastq.gz')]
    

#%% Map reads against IS6110
readparsing.mapreads(
    reads, target, 'reads_vs_IS', temp_dir, 'paf', args.cpus, k=9, m=10)


#%% Create read dictionary
from lib import classi_e_funzioni as cef
read_dict = cef.parse_paf(f'{temp_dir}/reads_vs_IS.paf')

# Traverse fastq and get read sequences. 
# Optionally write file with complete reads for reference mapping!
cef.add_seqs_to_read_dict(read_dict, reads, temp_dir)


#%% Cluster reads with cd-hit-est, separately for 5' and 3' side
anchor_clusters = {
    '5' : {},
    '3' : {}
}

for side in anchor_clusters:
    
    clusters = cef.cd_hit(
        f'{temp_dir}/anchors.{side}.fasta', 
        f'{temp_dir}/cd_hit_{side}')
    
    for cluster_id in clusters:
        
        anchor_cluster = cef.AnchorCluster(cluster_id, side)
        
        for read in clusters[cluster_id]:
            anchor_cluster.add_read(read, read_dict)  ### ambiguous read_dict!
            
        anchor_cluster.align_anchor_reads(temp_dir, args)
        anchor_cluster.get_cluster_consensus(temp_dir)
        anchor_clusters[side][cluster_id] = anchor_cluster
        
        
#%% Summarize output 
for side in anchor_clusters:
    for cluster_id in anchor_clusters[side]:
        print(cluster_id, len(anchor_clusters[side][cluster_id].reads))
        


#%% Optional: identify reference positions

# Map partially mapping reads against reference
readparsing.mapreads(
    [f'{temp_dir}/partially_mapping.fastq.gz'], reference, 
    'reads_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10)


"""
# Get split reads from reference alignment
splitreads = readparsing.getsplitreads(
    f'{temp_dir}/reads_vs_ref.bam', temp_dir, args.min_split_len, args.mapq)

# Extract partially mapping reads    
partially_mapping = readparsing.get_partially_mapping(f'{temp_dir}/reads_vs_IS.paf')

# Write to fastq
readparsing.subset_fastq(partially_mapping, reads, temp_dir)
anchors = readparsing.write_anchor_sequences(f"{temp_dir}/partially_mapping.fastq.gz", partially_mapping, temp_dir)
"""


#%% Create output



    
        