#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import os
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
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
        """
        Initialize a Read object with a read ID.
        
        Parameters
        ----------
        read_id : str
            ID of the read
            
        """
        self.read_id = read_id
        self.paf_rows = []  # To investigate secondary mappings
        
        self.anchor_coordinates = []
        self.target_coordinates = []
        
        self.anchor_seq = []
        self.target_seq = []
    
    def add_coordinates(self, paf_row, boundary_margin):
        """
        Add the coordinates of the anchor and the IS to the object, given a paf row.
        Use lists, since there might be multiple mappings
        
        Parameters
        ----------
        paf_row : list
            A list of strings, where each element is a field from a paf row.
            
        """
        self.paf_rows.append(paf_row)
        
        # Anchor part (read part that does not map against IS)
        query_start = int(paf_row[2])
        query_end = int(paf_row[3])
        query_len = int(paf_row[1])
        
        not_mapping = set(range(query_len)) - set(range(query_start, query_end))

        anchor_start = min(list(not_mapping))
        anchor_end = max(list(not_mapping))
        self.anchor_coordinates.append((anchor_start, anchor_end))

        # IS part
        target_start = int(paf_row[7])
        target_end = int(paf_row[8])
        target_len = int(paf_row[6])
        
        # To which side of the IS does the read map? 
        # Assuming that the start and end of the element are the same in the reads 
        # as in the provided sequence           
        if target_start <= boundary_margin:
            self.side = '5'
        elif target_end >= target_len - boundary_margin:
            self.side = '3'
        else:
            print('Check IS boundary conditions:', self.read_id, target_start, target_end)
        
        # Part of the IS covered by the read
        self.target_coordinates.append((target_start, target_end))
        
        
    
    def add_sequences(self, read):
        
        """
        Add sequences of the anchor and target parts to the read.
        
        Parameters
        ----------
        read : SeqRecord
            The read from which to extract the sequences.
        """
        for i, anchor in enumerate(self.anchor_coordinates):
             
            target = self.target_coordinates[i]
            
            anchor_start, anchor_end = anchor
            target_start, target_end = target
                        
            self.anchor_seq.append(
                SeqRecord(
                    read.seq[anchor_start:anchor_end]    ,
                    id=f'{read.id}_{i}',
                    name = '',
                    description=f'{anchor_start}-{anchor_end}'
                )    
            )
            
            self.target_seq.append(
                SeqRecord(
                    read.seq[target_start:target_end],
                    id=f'{read.id}_{i}',
                    name = '',
                    description=f'{target_start}-{target_end}'
                    )
                )
        

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
    
    
    def add_read(self, read_id, read_dict):
        """
        Add a read to the AnchorCluster instance.
        
        Parameters
        ----------
        read_id : str
            ID of the read to add.
        read_dict : dict
            Dictionary with read IDs as keys and anchor sequences as values.
        """
        anchor_rec = read_dict[read_id].anchor
        #anchor_rec.description = f'{self.side}_{self.cluster_nr}'
        self.reads.append(anchor_rec)
    
    
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
        
    def summarize_IS_coordinates(self, read_d):
        """ 
        
        
        """
        pass
                
        

    def get_cluster_consensus(self, temp_dir):
        
        """
        Compute the consensus sequence for the cluster based on alignment.

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
        
        consensus = consensus.strip('-').upper()
        
        self.prop_sites_with_mismatches = round(n_sites_with_mismatches / self.aln_len, 2)

        self.consensus = SeqRecord(
            Seq(consensus),
            id=self.cluster_id,
            name = '',
            description=''
            )
        




#%% Funzioni

def mapreads(fastq, ref, outpref, outpath, outfmt, cpus=1, k=15, m=40, more_args = []):
    """ Map Illumina reads against a reference. 

    Parameters
    ----------
    fastq : list
        Paths to fastq file(s).
    ref : str
        Path to reference against which reads are mapped.
    outpref : str
        Prefix for the output file.
    outpath : str
        Path to output folder.    
    outfmt : str
        Output format, either paf or bam.

    cpus : int, optional
        Number of CPUs. The default is 1.
    k : int, optional
        k-mer size for indexing in minimap2. The default is 15.
    m : int, optional
        Minimal chaining score. The default is 40.
    more_args : list, optional
        Additional arguments for minimap2. The default is [].

    Returns
    -------
    A paf or bam file with mapped reads, located in tmp/


    """

    cpus = str(cpus)

    if outfmt == 'bam':
  
        minimap_cmd = [
            'minimap2', '-a', '-x', 'sr', '-Y','-t', cpus, '-k', str(k), '-m', str(m)
        ] + more_args + [ref]

    elif outfmt == 'paf':
  
        minimap_cmd = [
            'minimap2', '-x', 'sr', '-c', '-o', f'{outpath}/{outpref}.paf', '-Y','-t', cpus, '-k', str(k), '-m', str(m)
        ] + more_args + [ref]

    else:
        return 'Error: specify output format'

    minimap_cmd += fastq

    if outfmt == 'bam':  # Get sorted and indexed bam file

        minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        subprocess.check_output(
            ('samtools', 'view', '-@', cpus, '-hu', '-o', f'{outpath}/{outpref}.raw.bam'),
            stdin=minimap.stdout, stderr=subprocess.DEVNULL)

        minimap.wait()

        subprocess.check_call(
            ('samtools', 'sort', '-@', cpus, f'{outpath}/{outpref}.raw.bam', '-o', f'{outpath}/{outpref}.bam'), stderr=subprocess.DEVNULL
            )

        subprocess.check_call(
            ("samtools", "index", f'{outpath}/{outpref}.bam')
            )

        subprocess.check_output(
            ('rm', f'{outpath}/{outpref}.raw.bam')
            )

    elif outfmt == 'paf':

        subprocess.run(minimap_cmd, stderr=subprocess.DEVNULL)
        
        
def parse_paf(paf_file, min_anchor_len=20, min_hit_len=20, boundary_margin=5):
    
    """  Traverse IS alignment to extract coordinates of IS and anchor read parts.
 
    Output a dictionary with read IDs, containing info about both the anchor and the IS part. 
     
    Complications:
        - nested insertions
        - close-by insertions
    """
    
    def is_overlapping(paf_row, min_anchor_len, min_hit_len, boundary_margin):
        """ 
        Test if read reaches into IS. 
        """
        read_len = int(paf_row[1])
        aln_len = int(paf_row[10])
        
        # Alignment too short
        if aln_len < min_hit_len:
            return False
        
        # Anchor part too short or read entirely in IS
        if (read_len - aln_len) < min_anchor_len:
            return False
        
        # Read do not map into IS, with margin of n bp
        target_start = int(paf_row[7])
        target_end = int(paf_row[8])
        target_len = int(paf_row[6])
        
        if not ( (target_start < boundary_margin) or (target_end > (target_len - boundary_margin)) ):
            return False
        
        else:
            return True

    read_d = {}

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
            
            if is_overlapping(fields,min_anchor_len, min_hit_len, boundary_margin):
                
                if read_id not in read_d:
                    read_d[read_id] = Read(fields)    
                                
                read_d[read_id].add_coordinates(fields, boundary_margin)
    
    return read_d


def add_seqs_to_read_dict(read_dict, reads, temp_dir):
    
    """
    Add sequences from FASTQ files to the read dictionary and optionally write them to FASTA files.

    This function processes a list of FASTQ files, extracting sequences for
    reads present in the provided read dictionary. The sequences are added to
    each read's corresponding entry in the dictionary. Optionally, the function
    can write the sequences to separate FASTA files for each side ('5' and '3').

    Parameters
    ----------
    read_dict : dict
        A dictionary where keys are read IDs and values are Read objects.
    reads : list
        A list of paths to FASTQ files containing the reads.
    temp_dir : str
        The path to the temporary directory where FASTA files will be written.

    """

    fasta_out = {
        '5' : [],
        '3' : [] 
        }
    
    for fastq in reads:
    
        with gzip.open(fastq, "rt") as fastq_handle:
        
            for read in SeqIO.parse(fastq_handle, "fastq"):
                
                if read.id in read_dict:
                    read_dict[read.id].add_sequences(read)
             
                    side = read_dict[read.id].side
                    for anchor_seq in read_dict[read.id].anchor_seq:
                        fasta_out[side].append(anchor_seq)

    for side in fasta_out:              
        with open(f'{temp_dir}/anchors.{side}.fasta', 'w') as fasta_handle:
            SeqIO.write(fasta_out[side], fasta_handle, 'fasta')

    
                       



def cd_hit(fasta_path, output_path):
    
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
        '-c', '0.95',
        '-o', output_path,
        '-sc', '1'
        ]
        
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


def parse_clusters(read_d,temp_dir, args):
    
    """
    Parse clusters of anchor sequences and compute their consensus.

    This function processes anchor sequences for both 5' and 3' sides by
    clustering them using `cd-hit-est`. For each cluster, it creates an
    `AnchorCluster` instance, adds reads to it, aligns the reads using MAFFT,
    and computes the consensus sequence.

    Parameters
    ----------
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

    anchor_clusters = {'5':{}, '3':{}}

    for side in anchor_clusters:
        
        clusters = cd_hit(
            f'{temp_dir}/anchors.{side}.fasta', 
            f'{temp_dir}/cd_hit_{side}')
        
        for cluster_id in clusters:
            
            anchor_cluster = AnchorCluster(cluster_id, side)
            
            for read in clusters[cluster_id]:  # read_id changed!
                read_id = read[:-2]
                read_index = int(read[-1])
                anchor_rec = read_d[read_id].anchor_seq[read_index]
                anchor_cluster.append(anchor_rec)
                
            anchor_cluster.align_anchor_reads(temp_dir, args)
            anchor_cluster.get_cluster_consensus(temp_dir)
            anchor_clusters[side][cluster_id] = anchor_cluster
        
    return anchor_clusters