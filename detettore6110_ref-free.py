#!/usr/bin/env python3

# %%

import os
import shutil
import tempfile

from Bio import SeqIO
from lib import classi_e_funzioni as cef
from lib import readparsing

class args:
    def __init__(self):
        self.reads = ['testing/reads.orygis.fastq.gz']
        self.target = 'resources/is_targets/IS6110.fasta'
        self.cpus = '4'
        self.reference = 'resources/reference/MTBC0_v1.1.fasta'

args = args()

working_dir = os.getcwd()    
reads = [os.path.abspath(x) for x in args.reads]
target = os.path.abspath(args.target)


#temp_dir = tempfile.mkdtemp()
# Keep files for developing...
temp_dir = 'testing/tmp'
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
os.mkdir(temp_dir)

# Convert input bam/cram to fastq
read_suffix = set([x.split('.')[-1] for x in reads]).pop()

if len(reads) == 1 and read_suffix in ['bam', 'cram', 'sam']:
    bamfile = reads[0]        
    reads = [readparsing.bam_to_fastq(bamfile, f'{temp_dir}/reads.fastq.gz')]

# Map reads against IS target
cef.mapreads(
    reads, target, 'reads_vs_IS', temp_dir, 'paf', args.cpus, k=9, m=10)


# Create read dictionary
read_d = cef.parse_paf(f'{temp_dir}/reads_vs_IS.paf')

# Traverse fastq and add anchor and hit parts of the reads
cef.add_seqs_to_read_dict(read_d, reads, temp_dir)


# Cluster reads with cd-hit-est and summarize clusters
clusters = cef.parse_clusters(read_d, temp_dir, args)   



#%% Optional: identify reference positions

readparsing.mapreads(
    [f'{temp_dir}/anchor_consensi.fasta'], args.reference, 
    'reads_vs_ref', temp_dir, 'bam', args.cpus, k=9, m=10)




        
#%% Summarize output 
for side in anchor_clusters:
    for cluster_id in anchor_clusters[side]:
        print(cluster_id, len(anchor_clusters[side][cluster_id].reads))












#%%
import pysam
bam = f'{temp_dir}/reads_vs_ref.bam'
pybam = pysam.AlignmentFile(bam, "rb")

for read in pybam.fetch():
    print(read.query_name, 
          read.reference_start, read.reference_end, 
          read.cigarstring, read.mapping_quality)
    

pybam.close()


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



    
        