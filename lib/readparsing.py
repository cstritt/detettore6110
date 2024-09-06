#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import gzip
import pysam
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
        
        cigar = read.cigartuples
        
        if len(cigar) > 1:
    
            first, last = cigar[0], cigar[-1]
            
            # what if both ends are clipped? like this, first should be called
            # more often
            if first[0] in (4,5):
                cliplen = first[1]

            elif last[0] in(4,5):
                cliplen = last[1]
                
            else:
                continue
            
            if cliplen < min_split_len:
                continue
            
            if read.mapping_quality < mapq:
                 continue

            # Write clipped to fasta
            name = read.query_name
            if read.is_read1:
                name = name +'/1'
            elif read.is_read2:
                name = name + '/2'
            
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
    # with open(f'{outpath}/partially_mapping.readIDs.txt', 'w') as f:
    #     for readid in partially_mapping:
    #         f.write(readid + '\n')
            
    # # Subset fastqs
    # for f in FASTQ:
    #     subprocess.run(
    #         ('seqtk', 'subseq', f, f'{outpath}/partially_mapping.readIDs.txt'),
    #         check=True
    #         )
    
    # Merge into one file
    
    fastq_out = []

    fasta_out = {
        '5' : [],
        '3' : []
        }
    
    anchor_d = {}
    
    for f in FASTQ: 
        
        with gzip.open(f, "rt") as fastq_handle:
        
            for read in SeqIO.parse(fastq_handle, "fastq"):
                
                if read.id in partially_mapping:
                    
                    fastq_out.append(read)
                    
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
             
    with gzip.open(f'{outpath}/partially_mapping.fastq.gz', 'wt') as fastq_handle:
        SeqIO.write(fastq_out, fastq_handle, 'fastq')
    
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

    samtools_cmd = ('samtools', 'fastq', '-0', outpath, bamfile)

    # Create a gzip process and pipe the output to it
    with open(outpath, "wb") as f:
        subprocess.run(samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout
        p = subprocess.Popen(["gzip"], stdin=subprocess.PIPE, stdout=f, stderr=subprocess.DEVNULL)
        p.communicate(subprocess.run(samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout)

    return outpath

