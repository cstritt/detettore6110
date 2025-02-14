#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pysam
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter
from sklearn.linear_model import LinearRegression

from collections import Counter


class AnchorCluster:

    def __init__(self,cluster_nr, side):
        """
        Initialize an AnchorCluster instance with a cluster number and side.

        Parameters
        ----------
        cluster_nr : intipython
        """

        self.cluster_nr = cluster_nr
        self.side = side
        self.cluster_id = f'{side}prime_{cluster_nr}'
        self.reads = []
        self.seqs = []
    
    def align_anchor_reads(self, temp_dir, args):
        
        """
        Align the reads in the AnchorCluster instance with MAFFT.
        
        Parameters
        ----------
        temp_dir : str
            Path to temporary directory.
        args : class
            Input arguments.
        """
        fasta_path = os.path.join(temp_dir, f'{self.cluster_id}.fasta')
        alignment_path = os.path.join(temp_dir, f'{self.cluster_id}.aligned.fasta')
        
        with open(fasta_path, 'w') as fasta_handle:
            SeqIO.write(self.seqs, fasta_handle, 'fasta')
        
        mafft_cmd = [
            'mafft',
            '--thread', str(args.cpus),
            '--adjustdirection',
            fasta_path
            ]
        
        subprocess.run(mafft_cmd, check=True, 
                    stdout=open(alignment_path, 'w'), 
                    stderr=subprocess.DEVNULL)
        
        os.remove(fasta_path)
        
        
    def summarize_IS_coordinates(self, read_d):
        """ Go through reads and store which positions of the IS are covered
        """
        
        self.target_cov = {}
        for read_id in self.reads:
            
            read = read_d[read_id]
            for coords in read.target_coordinates:
                start, end = coords
                for i in range(start, end+1):
                    if i not in self.target_cov:
                        self.target_cov[i] = 0
                    self.target_cov[i] += 1
                
        depth = [self.target_cov[i] for i in self.target_cov]
        
        pos = [[i] for i in self.target_cov]  # sort and use index rather than actual position
        
        self.target_lm = lm(depth, pos)
                

    def get_cluster_consensus(self, temp_dir):
        
        """
        Compute the consensus sequence for the cluster based on alignment.
        
        To assess the quality of the cluster, fit an lm(coverage~position). 
        Slope and intercept show if the coverage is increasing or decreasing as 
        expected, or if the cluster is messy and no trend in coverage is there.

        This method reads the alignment from a given file, calculates the
        consensus sequence, and determines the proportion of sites with
        mismatches. It also tracks the sequence depth at the start and end
        of the alignment.

        Parameters
        ----------
        alignment_path : str
            Path to the alignment file in FASTA format.

        Attributes
        ----------
        nr_reads : int
            Number of reads in the alignment.
        aln_len : int
            Length of the alignment.
        consensus : str
            Consensus sequence derived from the alignment.
        depth_start : int
            Depth of sequence coverage at the start of the alignment.
        depth_end : int
            Depth of sequence coverage at the end of the alignment.
        prop_sites_with_mismatches : float
            Proportion of sites with mismatches in the alignment.
        """
        
        alignment_path = os.path.join(temp_dir, f'{self.cluster_id}.aligned.fasta')

        aln = AlignIO.read(open(alignment_path), "fasta")
        aln_smry = AlignInfo.SummaryInfo(aln)
        
        self.nr_reads = len(aln)
        self.aln_len = aln.get_alignment_length()
        
        consensus = ''
        n_sites_with_mismatches = 0
        
        # For linear model
        depth = []
        position = []
        
        for i in range(self.aln_len):
            col = aln_smry.get_column(i)
            count_missing = col.count('-')
            count_present = self.nr_reads - count_missing
            
            depth.append(count_present)
            position.append([i+1])
            
            # Check on which side the "tail" of the alignment is
            if i == 0:
                self.depth_start = count_present
            if i == (self.aln_len-1):
                self.depth_end = count_present         
            
            # Get consensus base
            if count_present >= 3:
                counter = Counter(col.replace('-', ''))  # don't! Majority can be gap
                base = counter.most_common(1)[0][0]
                
                if len(counter) > 1:
                    n_sites_with_mismatches += 1
            else:
                base = '-'
                
            consensus += base
        
        consensus = consensus.replace('-', '').upper()
        
        self.prop_sites_with_mismatches = round(n_sites_with_mismatches / self.aln_len, 2)

        self.consensus = SeqRecord(
            Seq(consensus),
            id=self.cluster_id,
            name = '',
            description=''
            )
        
        self.anchor_lm = lm(depth, position)

    def add_reference_coordinates(self, ref_d):
        if self.cluster_id in ref_d:
            read = ref_d[self.cluster_id]
            self.ref = read.reference_name
            self.ref_start = read.reference_start
            self.ref_end = read.reference_end
            self.ref_strand = '-' if read.is_reverse else '+'
            self.ref_cigar = read.cigarstring
            self.ref_mapq = read.mapping_quality
        else:
            self.ref = 'NA'
            self.ref_start = 'NA'
            self.ref_end = 'NA'
            self.ref_strand = 'NA'
            self.ref_cigar = 'NA'
            self.ref_mapq = 'NA'
   
        
       
        
def cluster_anchors(read_d, l):
    """_summary_

    Args:
        read_d (_type_): _description_
        l (_type_): _description_
    """
    
    def group_strings_by_value(dictionary):
        """ Group a dictionary by values. Used here to cluster sequences by exact identity
        
        Args:
            dictionary (dict): keys are sequence IDs, values are sequences.

        Returns:
            dict: keys are sequences, values are sequence IDs
        """

        grouped = {}
        for key, values in dictionary.items():
            for value in values:
                grouped.setdefault(value, []).append(key)
        return grouped
    
    # Clustering, separate for 5' and 3'
    
    identifiers_5 = {}
    identifiers_3 = {}

    for read_id in read_d:
        for rec in read_d[read_id].anchor_seq:
            
            seq = str(rec.seq)
            seq_id = rec.id
            
            if read_d[read_id].side == '5':
                identifiers_5.setdefault(seq_id, []).append(seq[-l:])
            
            elif read_d[read_id].side == '3':
                identifiers_3.setdefault(seq_id, []).append(seq[:l])
                
    clusters_5 = group_strings_by_value(identifiers_5)
    clusters_3 = group_strings_by_value(identifiers_3)
    
    return {'5':clusters_5, '3':clusters_3}


def cd_hit(fasta_path, output_path, min_id = 0.99, ap=False):
    
    """
    Run cd-hit-est on a fasta file of anchor sequences and return a dictionary
    where keys are cluster numbers and values are lists of read IDs present in
    each cluster.

    Parameters
    ----------
    fasta_path : str
        The path to the fasta file to be clustered.
    output_path : str
        The path to the output file (without the .clstr extension).

    Returns
    -------
    clusters : dict
        A dictionary where keys are cluster numbers and values are lists of
        read IDs present in each cluster.
    """
    cd_hit = [
        'cd-hit-est',
        '-i', fasta_path,
        '-d', '0',
        '-c', str(min_id),
        '-o', output_path,
        '-sc', '1', 
        '-g', '1'
        ]
        
    if ap:
        cd_hit += ['-ap', '1']  # Doesn't seem to do anything...
        
    subprocess.run(cd_hit, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Parse output
    clusters = {}
    
    with open(f'{output_path}.clstr') as f:
        for line in f:
        
            if line.startswith('>'):
                cluster_nr = line.strip().split(' ')[-1]
                clusters[cluster_nr] = []
            else:
                read_id = line.strip().split(' ')[1][1:-3]
                clusters[cluster_nr].append(read_id)
    
    return clusters


def lm(target, features):
    """
    Perform linear regression on depth and positions data.

    This function fits a linear regression model to the given depth and positions
    data, and returns the intercept and slope of the fitted line.

    Parameters
    ----------
    depth : array-like
        The independent variable data (e.g., depth values).
    positions : array-like
        The dependent variable data (e.g., position values).

    Returns
    -------
    intercept : float
        The y-intercept of the regression line.
    slope : float
        The slope of the regression line.
    """

    m = LinearRegression()
    m.fit(features,target)
    intercept = round(float(m.intercept_),2)
    slope = round(float(m.coef_[0]),2)
    return intercept, slope
    

def parse_clusters(read_d, anchor_clusters, temp_dir, args):
    
    """
    Parse clusters of anchor sequences and compute their consensus.

    This function processes anchor sequences for both 5' and 3' sides by
    clustering them using `cd-hit-est`. For each cluster, it creates an
    `AnchorCluster` instance, adds reads to it, aligns the reads using MAFFT,
    and computes the consensus sequence.

    Parameters
    ----------
    read_d : dict
        A dictionary where keys are read IDs and values are Read objects.
    temp_dir : str
        Path to the temporary directory where intermediate files are stored.
    args : class
        Input arguments containing parameters such as number of CPUs for MAFFT.

    Returns
    -------
    anchor_clusters : dict
        A dictionary containing `AnchorCluster` objects for both 5' and 3' sides,
        keyed by cluster IDs.
    """

    cluster_d = {'5':{}, '3':{}}
    fasta_handle = open(f'{temp_dir}/anchor_consensi.fasta', 'w')

    for side in anchor_clusters:
        
        n = 0
        
        for seq in anchor_clusters[side]:
            
            # Skip short clusters
            if len(anchor_clusters[side][seq]) < args.min_cluster_size:
                continue
            
            ac = AnchorCluster(n, side)
            
            for read in anchor_clusters[side][seq]: 
                read_id = read[:-2]
                read_index = int(read[-1])
                anchor_rec = read_d[read_id].anchor_seq[read_index]
                    
                ac.seqs.append(anchor_rec)
                ac.reads.append(read_id)
                
            ac.align_anchor_reads(temp_dir, args)
            ac.get_cluster_consensus(temp_dir)
            ac.summarize_IS_coordinates(read_d)
            cluster_d[side][n] = ac
                
            if len(ac.consensus) > args.min_anchor_len:
                SeqIO.write(ac.consensus, fasta_handle, 'fasta')
        
            n += 1
            
    fasta_handle.close()
                
    return cluster_d


def find_overlaps(ref_aligned_anchors, tsd_len):
    
    """
    Find overlapping reads in a given sorted BAM file.

    Parameters
    ----------
    ref_aligned_anchors : str
        Path to the sorted BAM file containing anchor reads aligned to the reference.
    tsd_len : list
        List of possible target site duplication lengths.

    Returns
    -------
    overlaps : list
        List of tuples containing read IDs, start and end coordinates of overlapping reads,
        and the length of the overlap.
    """
    pybam = pysam.AlignmentFile(ref_aligned_anchors, "rb")
    read_ids = []
    coordinates = []

    for read in pybam.fetch():
        read_ids.append(read.query_name)
        coordinates.append((read.reference_start, read.reference_end))
    pybam.close()

    overlaps = []
    for i in range(len(coordinates)):
        for j in range(i + 1, len(coordinates)):
            start1, end1 = coordinates[i]
            start2, end2 = coordinates[j]
            overlap = end1 - start2
            if overlap in tsd_len:
                if read_ids[i][0] == read_ids[j][0]:  # make sure reads reach into opposite ends of the IS
                    continue
                overlaps.append((read_ids[i], end1, read_ids[j], start2, overlap))
                
    return overlaps


def add_ref_coordinates_to_clusters(cluster_d, ref_aligned_anchors):
    
    """
    Add reference coordinates to AnchorCluster objects from a sorted BAM file.

    Parameters
    ----------
    cluster_d : dict
        Dictionary containing AnchorCluster objects for both 5' and 3' sides,
        keyed by cluster IDs.
    ref_aligned_anchors : str
        Path to the sorted BAM file containing anchor reads aligned to the reference.

    Returns
    -------
    None

    """
    
    pybam = pysam.AlignmentFile(ref_aligned_anchors, "rb")
    
    ref_d = {}

    for read in pybam.fetch():
        ref_d[read.query_name] = read
    pybam.close()
    
    for side in cluster_d:
        for cluster_id in cluster_d[side]:
            cluster_d[side][cluster_id].add_reference_coordinates(ref_d)

