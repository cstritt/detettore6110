#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
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


def parse_paf(paf_file, min_anchor_len, min_hit_len):

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
    read_d = {}
    
    paf_header = [
        'query_name', 'query_length', 'query_start', 'query_end', 'strand',
        'target_name', 'target_length', 'target_start', 'target_end', 
        'residue_matches', 'block_length', 'mapping_quality'
    ]

    paf = pandas.read_csv(paf_file, sep='\t', usecols=range(12), header=None)
    paf.columns = paf_header
    
    for i, row in paf.iterrows():
        
        alignment_length = row['query_end'] - row['query_start']
        
        # Anchor not long enough or read entirely in IS
        if (row['query_length'] - alignment_length) < min_anchor_len:  
            continue
        
        # IS part not long enough
        if alignment_length < min_hit_len: 
            continue
        
        # Alignment does not beginn precisely at start or end of IS
        if not ((row['target_start'] == 0) or (row['target_end'] == row['target_length'])):
            continue
        
        # Query start or end is not in IS
        if not ((row['query_start'] == 0) or (row['query_end'] == row['query_length'])):
            continue

        if row['query_name'] not in read_d:
            read_d[row['query_name']] = Read(row['query_name'])
        
        read_d[row['query_name']].add_coordinates(row)

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



def get_partially_mapping(paf, min_anchor_len=20):
    """ 
    Find reads partially mapping against IS, given minimap2 paf output.    

    Parameters
    ----------
    paf : str
        Path to PAF output of minimap, containing reads that map to IS.
        
    min_anchor_len : int. Optional.
        Minimum length of the anchor/clipped read part. Default is 20.
        

    Returns
    -------
    
    partial : dict
    
    Dictionary with read IDs as keys. Values: the side of the IS to which the
    read maps, the start and the end coordinates of the anchor/clipped read part

    """
    
    partial = {}

    with open(paf) as f:
        
        for line in f:
            
            fields = line.strip().split('\t')
            
            query_len = int(fields[1])
            aln_len = int(fields[10])
            
            nonmap_len = query_len - aln_len
            
            if nonmap_len > min_anchor_len:
                
                readname = fields[0]
                
                # To which side of the IS does the read map?
                target_start = int(fields[7])
                target_end = int(fields[8])
                target_len = int(fields[6])
                
                query_start = int(fields[2])
                query_end = int(fields[3])
                
                if target_start == 0:
                    side = '5'

                elif target_end == target_len:
                    side = '3'
                    
                else:
                    continue
                
                if query_start == 0:
                    anchor_start = query_end
                    anchor_end = query_len
                
                elif query_end == query_len:
                    anchor_start = 0
                    anchor_end = query_start       
                    
                else:
                    continue
                
                partial[readname] = (side, anchor_start, anchor_end)
                        
    return partial

        
def getsplitreads(bam, outpath, min_split_len=15, mapq=1):
    """ Extract splitreads and write split parts to fasta file.
    

    Parameters
    ----------
    bam : str
        Path to bam file.
    outpath : str
        Path to output directory.
    min_len : int, optional
        Minimum length of the split part of the read. The default is 10.
    mapq : int, optional
        Minimum mapping quality of the split read. The default is 1.

    Returns
    -------
    splitreads : dict
        A dictionary with read IDs as keys and pysam read objects as values.
        
    Also writes a fasta file containing the clipped parts of reads at 
    outpath/softclipped.fasta.

    """

    pybam = pysam.AlignmentFile(bam, "rb")

    splitreads = {}
    
    splitseqs = []

    for read in pybam.fetch():
        
        if read.is_unmapped:
            continue
        
        cigar = read.cigartuples
        
        if len(cigar) == 2:
    
            first, last = cigar[0], cigar[-1]
            # what if both ends are clipped? like this, first should be called
            # more often
            if first[0] == 4:  # 4 stands for soft-clipped
                cliplen = first[1]

            elif last[0] == 4:
                cliplen = last[1]
                
            else:
                continue
            
            if (cliplen < min_split_len) or (read.mapping_quality < mapq):
                continue

            # Write clipped to fasta
            name = read.query_name
            if read.is_read1:
                name = name +'/1'
            elif read.is_read2:
                name = name + '/2'
                
            # Reset read query name: add suffix, which is sometimes removed by minimap2
            read.query_name = name
            
            seqrec = SeqRecord(
                Seq(clip_seq(read.query_sequence, cigar)), 
                id = name, 
                name="", 
                description="")
            
            if name not in splitreads:
                splitreads[name] = []
             
            splitreads[name].append(read)  
            
            splitseqs.append(seqrec)
        
    pybam.close()
    
    SeqIO.write(splitseqs, outpath + '/softclipped.fasta', 'fasta')

    return splitreads


def clip_seq(seq, cigar):
    """     Return clipped part of sequence, as indicated in cigar.
    Assumes that there is only one clipped part at the beginning or end of 
    the read. Cigar input parameter is a pysam cigartuples object.
    

    Parameters
    ----------
    seq : str or Seq
        The sequence including soft clipped bases.
    cigar : cigartuple
        The cigar of the sequence in pysams cigartuple class.

    Returns
    -------
    clipseq : Seq
        The softclipped part of the read.

    """
    
    if cigar[0][0] == 4:
        end = cigar[0][1]
        clipseq = seq[:end]
    elif cigar[-1][0] == 4:
        strt = sum([x[1] for x in cigar[:-1]])
        end = strt + cigar[-1][1]
        clipseq = seq[strt:end]
    return clipseq


def subset_fastq(partially_mapping, FASTQ, outpath):
    """
    Extract the partially mapping reads from the original fastq. 
    
    Modify: first subset with seqtk (much faster!), then create fasta with 
    anchor sequences
    
    
    Parameters
    ----------
    partially_mapping : dict
        Output of the get_partially_mapping function.
    FASTQ : list
        Original FASTQ file(s).
    args : class
        Input arguments.

    Returns
    -------
    
    anchor_d : dict
        A dictionary with read IDs as keys and anchor sequences as values
    
    One (SE) or two (PE) fastq files containing the IS-mapping reads,
    two fasta files containing the anchors of the 5' and the 3' sides, written
    to outpath.

    """
    
    # Write read IDs to file
    read_ids = f'{outpath}/partially_mapping.readIDs.txt'
    with open(read_ids, 'w') as f:
        for readid in partially_mapping:
            f.write(readid + '\n')
            
    # PE reads
    if len(FASTQ) == 2:
        
        # Subsample the forward reads
        forward_subset = f"{outpath}/subsampled_forward.fastq"
        subprocess.run([
            "seqtk", "subseq", FASTQ[0], read_ids
        ], check=True, stdout=open(forward_subset, "w"))

        # Subsample the reverse reads
        reverse_subset = f"{outpath}/subsampled_reverse.fastq"
        subprocess.run([
            "seqtk", "subseq", FASTQ[1], read_ids
        ], check=True, stdout=open(reverse_subset, "w"))

        # Merge the subsampled reads
        merged_fastq = f"{outpath}/partially_mapping.fastq"
        subprocess.run([
            "cat",  forward_subset, reverse_subset
        ], check=True, stdout=open(merged_fastq, "w"))
    
    # SE reads
    elif len(FASTQ) == 1:
        se_subset = f"{outpath}/partially_mapping.fastq"
        subprocess.run([
            "seqtk", "subseq", FASTQ[0], read_ids
        ], check=True, stdout=open(se_subset, "w"))
        
    # Compress
    subprocess.run([
        "gzip", "-f", f"{outpath}/partially_mapping.fastq"
        ], check=True, stdout=open(f"{outpath}/partially_mapping.fastq.gz", "wb"))
        
                
def write_anchor_sequences(FASTQ, partially_mapping, outpath):
    """
    Write the anchor sequences to fasta file.

    Parameters
    ----------
    FASTQ : str
        Path to fastq file. This is a single file resulting form subset_fastq above.
    partially_mapping : dict
        Dictionary with read IDs as keys and the side, start and end coordinates
        of the anchor as values.
    outpath : str
        Path to output folder.

    Returns
    -------
    anchor_d : dict
        Dictionary with read IDs as keys and anchor sequences as values.

    """
    fasta_out = {
        '5' : [],
        '3' : []
        }
    
    anchor_d = {}
    
    with gzip.open(FASTQ, "rt") as fastq_handle:
    
        for read in SeqIO.parse(fastq_handle, "fastq"):
            
            if read.id in partially_mapping:
                
                side = partially_mapping[read.id][0]
                
                anchor_start = partially_mapping[read.id][1]
                anchor_end = partially_mapping[read.id][2]
                anchor_seq = read.seq[anchor_start:anchor_end]
                
                fasta_rec = SeqRecord(
                    anchor_seq,
                    id=read.id,
                    name = read.id,
                    description=side + '_' + str(anchor_start) + '-' + str(anchor_end)
                    )

                fasta_out[side].append(fasta_rec)
                anchor_d[read.id] = anchor_seq
                
    for side in fasta_out:
        with open(f'{outpath}/anchors.{side}.fasta', 'w') as fasta_handle:
            SeqIO.write(fasta_out[side], fasta_handle, 'fasta')

    return anchor_d


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

