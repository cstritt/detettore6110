#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Module to estimate IS copy numbers from split read clusters, without
mapping to reference

Approach:
    - cluster anchor read parts with cd-hit-est, separately for 5' and 3' side
    - CN = min(len(n_clusters_left), len(n_clusters_right))
    
"""

import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cluster_anchors(fasta, anchors, outpath):
    """
    Use cd-hit to cluster anchor reads, then assemble the reads of a cluster
    into a single sequence with CAP3.
    
    cd-hit parameters don't seem to greatly affect the results (word size, 
    minimum identity ...)
    

    Parameters
    ----------

    fasta : list
        List of two fasta files, for the 5' and the 3' anchors
    anchors : dict
        Dictionary with read ID as key and sequence as value, as output by
        subset_fastq.


    Returns
    -------    
    
    A fasta file with one assembled consensus sequence for each anchor cluster.
    
    """
    
    clusters = {
        '5' : {},
        '3' : {}
        }
    
    for f in fasta:
        
        side = f.split('.')[-2]
        
        cd_hit = [
            'cd-hit-est',
            '-i', f,
            '-d', '0',
            '-c', '0.95',
            '-o', f'{outpath}/cd_hit_{side}',
            '-sc', '1']
        
        subprocess.run(cd_hit, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # Collect anchor sequences for each cluster
        with open(f'{outpath}/cd_hit_{side}.clstr') as f:
            for line in f:
                
                if line.startswith('>'):
                    
                    cluster_nr = line.strip().split(' ')[-1]
                    clusters[side][cluster_nr] = []
                    
                else:
                    read_id = line.strip().split(' ')[1][1:-3]
                    anchor_seq = anchors[read_id]
                    anchor_rec = SeqRecord(
                        anchor_seq, 
                        id=read_id,
                        name='',
                        description = side +'_' + cluster_nr
                        )
                    
                    clusters[side][cluster_nr].append(anchor_rec)
                    
    return clusters



def get_copy_number(clusters, min_cluster_size=10):
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


def assemble_cluster_consensi(clusters, params):
    """ Use cap3 to assemble the reads of a cluster into a consensus sequence.
    
    Cave: assembled sequences are often few bp short at the clipped end! So
    better map the original reads.
    
    Better: use mafft to align!

    
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
            cap3_cmd = ['/home/cristobal/programs/CAP3/cap3', params.tmp + '/cluster.tmp.fasta']
            
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


