#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Module to estimate IS copy numbers from split read clusters, without
mapping to reference

Approach:
    - cluster anchor read parts with cd-hit-est, separately for 5' and 3' side
    - CN = min(len(n_clusters_left), len(n_clusters_right))
    
"""

import subprocess
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter


class AnchorCluster:
    
    """
    Info to add: to which part of the IS the reads map...
    
    
    cluster_id, side, nr_reads, depth_start, depth_end, 
    prop_sites_with_mismatches, len(consensus), consensus
    
    """
    def __init__(self,cluster_nr, side):
        self.cluster_nr = cluster_nr
        self.side = side
        self.cluster_id = f'{side}prime_{cluster_nr}'
        self.reads = []
    
    
    def add_read(self, read_id, read_dict):
        anchor_seq = read_dict[read_id]
        
        anchor_rec = SeqRecord(
            anchor_seq, 
            id=read_id,
            name='',
            description = f'{self.side}_{self.cluster_nr}'
            )
        
        self.reads.append(anchor_rec)
    
    
    def align_anchor_reads(self, temp_dir, args):
        
        fasta_path = os.path.join(temp_dir, f'{self.cluster_id}.fasta')
        alignment_path = os.path.join(temp_dir, f'{self.cluster_id}.fasta')
        
        with open(fasta_path, 'w') as fasta_handle:
            SeqIO.write(self.reads, fasta_handle, 'fasta')
        
        mafft_cmd = [
            'mafft',
            '--thread', args.cpus,
            '--adjustdirection',
            fasta_path
            ]
        
        subprocess.run(mafft_cmd, check=True, 
                    stdout=open(alignment_path, 'w'), 
                    stderr=subprocess.DEVNULL)
        
        
    def get_cluster_consensus(self, alignment_path):
        
        aln = AlignIO.read(open(alignment_path), "fasta")
        aln_smry = AlignInfo.SummaryInfo(aln)
        
        self.nr_reads = len(aln)
        self.aln_len = aln.get_alignment_length()
        
        consensus = ''
        n_sites_with_mismatches = 0
        
        for i in range(self.aln_len):
            col = aln_smry.get_column(i)
            count_missing = col.count('-')
            count_present = self.nr_reads - count_missing
            
            # Check on which side the "tail" of the alignment is
            if i == 0:
                self.depth_start = count_present
            if i == (self.aln_len-1):
                self.epth_end = count_present         
            
            # Get consensus base
            if count_present >= 3:
                counter = Counter(col.replace('-', ''))
                base = counter.most_common(1)[0][0]
                if len(counter) > 1:
                    n_sites_with_mismatches += 1
            else:
                base = '-'
                
            consensus += base
        
        self.consensus = consensus.strip('-').upper()
        self.prop_sites_with_mismatches = round(n_sites_with_mismatches / self.aln_len, 2)


def cd_hit(fasta_path, output_path):
    
    cd_hit = [
        'cd-hit-est',
        '-i', fasta_path,
        '-d', '0',
        '-c', '0.95',
        '-o', output_path,
        '-sc', '1'
        ]
        
    subprocess.run(cd_hit, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Parse output
    clusters = {}
    
    with open(output_path) as f:
        for line in f:
        
            if line.startswith('>'):
                cluster_nr = line.strip().split(' ')[-1]
                clusters[cluster_nr] = []
            else:
                read_id = line.strip().split(' ')[1][1:-3]
                clusters[cluster_nr].append(read_id)
    
    return clusters


def infer_copy_number(clusters, min_cluster_size=10):
    """

    Parameters
    ----------
    clusters : dict
        Output of cluster_anchors above.
    min_cluster_size : int
        Do not consider clusters with fewer than x reads.
    params : class
        Detettore parameters (mise_en_place class).

    Returns
    -------
    copy_number : int
        IS copy number.

    """

    n_clusters = {
        }

    for side in clusters:
        n_clusters[side] = 0
        for cl in clusters[side]:
            if len(clusters[side][cl]) < min_cluster_size:
                continue
            n_clusters[side] += 1

    copy_number = min([n_clusters['5'], n_clusters['3']])
    return copy_number   



def main(path_to_5prime_anchors, path_to_3prime_anchors, anchor_dict, temp_dir, args):

    clusters = {
        '5' :[],
        '3' : []
        }
    
    for side, f in map(
        ['5', '3'],
        [path_to_5prime_anchors, path_to_3prime_anchors]):
        
        # Cluster reads with cd-hit
        cd_hit_out = f'{temp_dir}/cd_hit_{side}'
        clusters = cd_hit(f, cd_hit_out)
        
        for cluster_nr in clusters:
            anchor_cluster = AnchorCluster(cluster_nr, side)
            for read in clusters[cluster_nr]:
                anchor_cluster.add_read(read, anchor_dict)

            anchor_cluster.align_anchor_reads(temp_dir, args)
            anchor_cluster.get_cluster_consensus()




#%% NOT USED

def assemble_cluster_consensi(clusters, params):
    """ Use cap3 to assemble the reads of a cluster into a consensus sequence.
    
    NOT USED, the mafft approach offers more control.
    
    Cave: assembled sequences are often few bp short at the clipped end! So
    better map the original reads.
    
    OR: use mafft to align!

    
    Parameters
    ----------
    clusters : dict
    
    
    Returns
    -------
    None.

    """
        
    # # Assemble the reads in each cluster
    cluster_consensi = {}
    
    for side in clusters:
        for cluster_nr in clusters[side]:
            
            cluster_id = side + '_' + cluster_nr
            
            with open(params.tmp + '/cluster.tmp.fasta', 'w') as fasta_handle:
                SeqIO.write(clusters[side][cluster_nr], fasta_handle, 'fasta')
                
            # Assemble
            cap3_cmd = ['cap3', params.tmp + '/cluster.tmp.fasta']
            
            subprocess.run(cap3_cmd, check=True, stderr=subprocess.DEVNULL)
            
            contigs = [seq_record.seq for seq_record in SeqIO.parse(params.tmp + '/cluster.tmp.fasta.cap.contigs', "fasta")]
            
            cluster_consensi[cluster_id] = contigs
            
        
    with open(params.tmp + '/cluster_consensi.fasta', 'w') as f:
           
        for cluster_id in cluster_consensi:
            
            contigs = cluster_consensi[cluster_id]
            
            # No consensus sequence
            if len(contigs) == 0:
                continue
            
            rec = SeqRecord(contigs[0], id = cluster_id, name='', description='')
        
            SeqIO.write(rec, f, 'fasta')
        
    return cluster_consensi
