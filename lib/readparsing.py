#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import os
import pandas
import pysam
import subprocess
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
            
    def add_coordinates(self, paf_row):
        """
        Add the coordinates of the anchor and the IS to the object, given a paf row.
        Use lists, since there might be multiple mappings
        
        In PAF format, the coordinates are 0-based, bed-like, 
        start closed and end open, meaning that the end coordinate is 
        not part of the alignment! The same coordinate format is used
        below and in the add_sequences function.
                
        Parameters
        ----------
        paf_row : list
            A list of strings, where each element is a field from a paf row.
            
        """
        self.paf_rows.append(paf_row.tolist())
        
        # Anchor part (read part that does not map against IS)
        self.query_strand = paf_row['strand']
        
        if paf_row['query_start'] == 0:
            anchor_start = paf_row['query_end']
            anchor_end = paf_row['query_length']
            
        elif paf_row['query_end'] == paf_row['query_length']:
            anchor_start = 0
            anchor_end = paf_row['query_start']
            
        else:
            warnings.warn(f'Check anchor coordinates for {paf_row['query_name']}\n')
            
        self.anchor_coordinates.append((anchor_start, anchor_end))

        # IS part
        # To which side of the IS does the read map? 
        # Assuming that the start and end of the element are the same in the reads 
        # as in the provided sequence           
        if paf_row['target_start'] == 0:
            self.side = '5'
        elif paf_row['target_end'] == paf_row['target_length']:
            self.side = '3'
        else:
            warnings.warn(f'Check IS boundaries for {paf_row["query_name"].tostring()}')
        
        # Part of the IS covered by the read
        self.target_coordinates.append((paf_row['target_start'], paf_row['target_end']))
        

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
            
            # Reorient reads such that they all begin with the IS overlapping part
            # This will allow more stric clustering with cd-hiz-est (-ap)
            anchor_seq = read.seq[anchor_start:anchor_end]
            read_id = f'{read.id}_{i}'
            
            if self.query_strand == '-':
                anchor_seq = anchor_seq.reverse_complement()                
                        
            self.anchor_seq.append(
                SeqRecord(
                    anchor_seq,
                    id=read_id,
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
            
                
class Reads:
    def __init__(self, args, temp_dir):
        
        # Input: convert to fastq if it's bam/cram/sam
        self.fastq = [os.path.abspath(x) for x in args.reads]
        read_suffix = set([x.split('.')[-1] for x in self.fastq]).pop()

        if len(self.fastq) == 1 and read_suffix in ['bam', 'cram', 'sam']:
            bamfile = self.fastq[0]        
            self.fastq = [bam_to_fastq(bamfile, f'{temp_dir}/reads.fastq.gz')]
            
        self.min_anchor_length = args.min_anchor_len
        self.min_hit_length = args.min_hit_len

        # Output 
        self.read_d = {}
        
        

    def parse_paf(self, paf_file, min_anchor_length, min_hit_length):        
        """
        Parse a minimap2 PAF file and return a dictionary of read IDs as keys 
        and Read objects as values.
        
        Parameters
        ----------
        paf_file : str
            Path to PAF output of minimap2.
            
        min_anchor_len : int
            Minimum length of the anchor (i.e. the non-aligned read part).
            
        min_hit_len : int
            Minimum length of the hit (i.e. the aligned part).
        
        Returns
        -------
        read_d : dict
            Dictionary with read IDs as keys and Read objects as values.
        """
        
        paf_header = [
            'query_name', 'query_length', 'query_start', 'query_end', 'strand',
            'target_name', 'target_length', 'target_start', 'target_end', 
            'residue_matches', 'block_length', 'mapping_quality'
        ]

        # Make sure file is not empty
        if os.stat(paf_file).st_size != 0:
             
            paf = pandas.read_csv(paf_file, sep='\t', usecols=range(12), header=None)
            paf.columns = paf_header
            
            for i, row in paf.iterrows():
                
                alignment_length = row['query_end'] - row['query_start']
                
                # Anchor not long enough or read entirely in IS
                if (row['query_length'] - alignment_length) < min_anchor_length:  
                    continue
                
                # IS part not long enough
                if alignment_length < min_hit_length: 
                    continue
                
                # Alignment does not beginn precisely at start or end of IS
                if not ((row['target_start'] == 0) or (row['target_end'] == row['target_length'])):
                    continue
                
                # Query start or end is not in IS
                if not ((row['query_start'] == 0) or (row['query_end'] == row['query_length'])):
                    continue

                if row['query_name'] not in self.read_d:
                    self.read_d[row['query_name']] = Read(row['query_name'])
                
                self.read_d[row['query_name']].add_coordinates(row)
    
    
    def add_seqs_to_read_dict(self, reads, temp_dir):
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
                    
                    if read.id in self.read_d:
                        self.read_d[read.id].add_sequences(read)
                
                        side = self.read_d[read.id].side
                        for anchor_seq in self.read_d[read.id].anchor_seq:
                            fasta_out[side].append(anchor_seq)

        for side in fasta_out:              
            with open(f'{temp_dir}/anchors.{side}.fasta', 'w') as fasta_handle:
                SeqIO.write(fasta_out[side], fasta_handle, 'fasta')


def mapreads(fastq, ref, outpref, outpath, outfmt, cpus=1, k=15, m=40):
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

    Returns
    -------
    A paf or bam file with mapped reads, located in tmp/


    """

    cpus = str(cpus)

    if outfmt == 'bam':
  
        minimap_cmd = [
            'minimap2', '-a', '-x', 'sr', '-Y','-t', cpus, '-k', str(k), '-m', str(m), ref
        ]

    elif outfmt == 'paf':
  
        minimap_cmd = [
            'minimap2', '-x', 'sr', '-c', '-o', f'{outpath}/{outpref}.paf', '-Y','-t', cpus, '-k', str(k), '-m', str(m), ref
        ]

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


def bam_to_fastq(bamfile, outpath):
    """
    Convert bam/cram to fastq using samtools

    -o option requires samtools 1.11!

    Parameters
    ----------
    bamfile : str
        Path to bam/cram file.

    Returns
    -------
    Writes fastq to outpath.

    """

    samtools_cmd = ('samtools', 'fastq', bamfile)

    # Create a gzip process and pipe the output to it
    with open(outpath, "wb") as f:
        subprocess.run(samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout
        p = subprocess.Popen(["gzip"], stdin=subprocess.PIPE, stdout=f, stderr=subprocess.DEVNULL)
        p.communicate(subprocess.run(samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout)

    return outpath

