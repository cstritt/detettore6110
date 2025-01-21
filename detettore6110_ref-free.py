#%%
import os
import tempfile

from lib import classi_e_funzioni as cef
from lib import readparsing

#%% mise en place
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


#%% Convert input bam/cram to fastq

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
cef.add_seqs_to_read_dict(read_dict, reads, temp_dir)











#%% Cluster reads with cd-hit-est, separately for 5' and 3' side
anchor_clusters = {
    '5' : {},
    '3' : {}
}

for side in anchor_clusters:
    clusters = cef.cd_hit(f'{temp_dir}/anchors.{side}.fasta', f'{temp_dir}/cd_hit_{side}', side)
    for cluster_id in clusters:
        anchor_cluster = cef.AnchorCluster(cluster_id, side)
        for read in clusters[cluster_id]:
            anchor_cluster.add_read(read, read_dict)
        anchor_cluster.align_anchor_reads(temp_dir, args)
        anchor_cluster.get_cluster_consensus()
        
    

#%% Create output file



    
        
        
#%% Assemble the reads in each cluster: mafft approach

""" 
Aim: create summary table with 

cluster_id, side, nr_reads, consensus 

Modify get_partially_mapping, so it also returns info about the IS:
    - which parts of the IS, how much
    - sequence identity
    
Output:
    - cluster_id
    - side
    - nr_reads
    - depth_start
    - depth_end
    - prop_sites_with_mismatches
    - len(consensus)
    - consensus
    

"""



import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter

min_depth = 3

cluster_consensi = {}

for side in anchor_clusters:
    for cluster_nr in anchor_clusters[side]:
        
        cluster_id = f'{side}prime_{cluster_nr}'
        
        with open(os.path.join(temp_dir, f'{cluster_id}.fasta'), 'w') as fasta_handle:
            SeqIO.write(anchor_clusters[side][cluster_nr], fasta_handle, 'fasta')
         
        mafft_cmd = [
            'mafft',
            '--thread', args.cpus,
            '--adjustdirection',
            os.path.join(temp_dir, f'{cluster_id}.fasta')
            ]
        
        subprocess.run(mafft_cmd, check=True, 
                       stdout=open(os.path.join(temp_dir, f'{cluster_id}.aligned.fasta'), 'w'), 
                       stderr=subprocess.DEVNULL)
        

        # Get consensus from alignment
        aln = AlignIO.read(open(os.path.join(temp_dir, f'{cluster_id}.aligned.fasta')), "fasta")
        nr_reads = len(aln)
        aln_len = aln.get_alignment_length()
        
        # Create consensus
        aln_smry = AlignInfo.SummaryInfo(aln)
        consensus = ''
        n_sites_with_mismatches = 0
        
        for i in range(aln_len):
            col = aln_smry.get_column(i)
            count_missing = col.count('-')
            count_present = nr_reads - count_missing
            prop_missing = count_missing / nr_reads
            
            # To check on which side the "tail" of the alignment is
            if i == 0:
                depth_start = count_present
            if i == (aln_len-1):
                depth_end = count_present         
            
            # Get consensus base            
            if count_present >= 3:
                counter = Counter(col.replace('-', ''))
                base = counter.most_common(1)[0][0]
                if len(counter) > 1:
                    n_sites_with_mismatches += 1
            else:
                base = '-'
                
            consensus += base                
        consensus = consensus.strip('-').upper()
        prop_sites_with_mismatches = round(n_sites_with_mismatches / aln_len, 2)
        
        row = [cluster_id, side, nr_reads, depth_start, depth_end, prop_sites_with_mismatches, len(consensus), consensus]    
        print(row)
        
        os.remove(os.path.join(temp_dir, f'{cluster_id}.fasta'))
