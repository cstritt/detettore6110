#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import os
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter


#%% Two main classes

class Read:
    """ 
    Class for reads that reach into the target IS. 
    Storing information for both the anchor part and the
    part mapping against the IS. 
    
    """
    def __init__(self, read_id):
        self.read_id = read_id
    
    def add_coordinates(self, paf_row):
        """
    
        query = read
        target = IS sequence
        
        """
        # Anchor part (read part that does not map against IS)
        query_start = int(paf_row[2])
        query_end = int(paf_row[3])
        query_len = int(paf_row[1])
        
        not_mapping = set(range(query_len)) - set(range(query_start, query_end))
        self.anchor_start = min(list(not_mapping))
        self.anchor_end = max(list(not_mapping))
        
            
        # IS part
        self.target_start = int(paf_row[7])
        self.target_end = int(paf_row[8])
        target_len = int(paf_row[6])
        
        # To which side of the IS does the read map? 
        # Assuming that the start and end of the element are the same in the reads 
        # as in the provided sequence           
        if self.target_start == 0:
            self.side = '5'
        elif self.target_end == target_len:
            self.side = '3'
        else:
            self.side = 'NA'  # What happened here?
            
        # Part of the IS covered by the read
        self.is_range = set(range(self.target_start, self.target_end))
        
    
    def add_sequences(self, read):
        
        anchor_seq = read.seq[self.anchor_start:self.anchor_end]        
        target_seq = read.seq[self.target_start:self.target_end]
                
        self.anchor = SeqRecord(
            anchor_seq,
            id=read.id,
            name = read.id,
            description=f'anchor_{self.side}_{self.anchor_start}-{self.anchor_end}'
            )
        
        self.targetpart = SeqRecord(
            target_seq,
            id=read.id,
            name = read.id,
            description=f'target_{self.side}_{self.target_start}-{self.target_end}' 
        )
        

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
    
    
class Summary:
    def __init__(self):
        self.reads = []
        self.clusters = []  





#%% Funzioni
def parse_paf(paf_file, min_anchor_len=20, min_hit_len=20):
    
    """  Traverse IS alignment to extract coordinates of IS and anchor read parts.
 
    Output a dictionary with read IDs, containing info about both the anchor and the IS part. 
     
    Complications:
        - nested insertions
        - close-by insertions
         
        
    """
    
    read_dict = {}

    with open(paf_file) as f:
        
        for line in f:
            
            fields = line.strip().split('\t')
            query_len = int(fields[1])
            aln_len = int(fields[10])
            
            if (query_len - aln_len) < min_anchor_len:  # anchor not long enough
                continue
            
            if aln_len < min_hit_len:  # IS part not long enought
                continue
            
            read_id = fields[0]
            
            if read_id in read_dict:
                print(f'Warning: {read_id} already in read_dict')
                continue
            
            read = Read(read_id)
            read.add_coordinates(fields)
            read_dict[read_id] = read
    
    return read_dict


def add_seqs_to_read_dict(read_dict, reads, temp_dir, write_fasta=True):
    
    fasta_out = {
        '5' : [],
        '3' : [], 
        'NA' : []  
        }
    
    for fastq in reads:
        
        print(fastq)
    
        with gzip.open(fastq, "rt") as fastq_handle:
        
            for read in SeqIO.parse(fastq_handle, "fastq"):
                
                if read.id in read_dict:
                    read_dict[read.id].add_sequences(read)
             
                    side = read_dict[read.id].side
                    fasta_out[side].append(read_dict[read.id].anchor)

    if write_fasta:
        for side in ['5', '3']:  
            with open(f'{temp_dir}/anchors.{side}.fasta', 'w') as fasta_handle:
                SeqIO.write(fasta_out[side], fasta_handle, 'fasta')

        
                       



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
