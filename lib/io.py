#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Module for writing output

"""

import bisect
import sys
import numpy
import os
import pandas
import pysam
import re
import shutil

from Bio import SeqIO
from collections import Counter
    
    
def exit_handler(args, temp_dir):
    """ Cleanup after program finish. If --keep is given, copy contents of 
    temporary directory to working directory before deleting it"""
    if args.keep:  # Copy contents of temporary to output directory
        shutil.copytree(temp_dir, os.path.join(args.outpath, args.prefix + '_intermediate_files'))
        
    shutil.rmtree(temp_dir)
    

def gene_overlap(position, annotation, chromosome_length):
    """  Given a genomic position, return the genomic context as given by 
    a gff annotation. Uses regex to extract strings after gene= and locus_tag=
    
    Issue:
        - in the h37rv annotation features can be two indices apart
        - for some features no regex is found
    
    
    Parameters
    ----------
    positions : DataFrame
        Pandas data frame containing the gff annotation.
    annotation : str
        Path to annotation file in gff format.
    chromosome_length : int
        Length of the reference genome. Used to get distance to dnaA when insertion is at the very end. 

    Returns
    -------
    A list with one element for each position, containing the gene name 
    and the distance to the gene. If the position is in a gene, this is one
    name and distance 0; if the position is between genes, this is two genes 
    and two distances separated by ;   

    """
    
    def gene_id_regex(patterns, gff_attributes):
        """ Extract gene name and locus tag (or any pattern) from the
        atrribute column of a gff
        
        Parameters
        ----------
        patterns : list
            List containing regex patterns to extract.
        gff_attributes : str
            Single gff attribute entry.

        Returns
        -------
        out : str
            Matches separated by a comma.

        """
        
        for pattern in patterns:
            match = re.search(pattern, gff_attributes)
            if match:
                return match.group(1)
    
    
    regex_patterns = [
        r";gene=([^;\n]+)", 
        r"locus_tag=([^;\n]+)", 
        r"mobile_element_type=([^;\n]+)"
        ]
    
    # Find the closest intervals on either side of each position
    idx_s = bisect.bisect_right(annotation['start'], position)
    idx_e = bisect.bisect_right(annotation['end'], position)
    
    # Overlapping feature
    if (idx_e == idx_s - 1) or (idx_e == idx_s -2):
                
        gene_info = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_e])
        
        return [gene_info, '0']
                                
    # Inbetween features
    elif idx_s == idx_e:
        
        gene_info_5 = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_s - 1])
        
        dist_to_5 = position - annotation['end'][idx_s - 1]
        
        # Insertion after last gene
        if idx_s == len(annotation):
            idx_s = 0
        
        gene_info_3 = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_s])
        
        if idx_s == 0:
            dist_to_3 = chromosome_length - position
        else:
            dist_to_3 = annotation['start'][idx_s] - position
        
        return [f'{gene_info_5};{gene_info_3}', f'{dist_to_5};{dist_to_3}']


def write_cluster_output(cluster_d, args):
    
    """
    Write output for anchor clusters to stdout.

    Parameters
    ----------
    cluster_d : dict
        Dictionary containing AnchorCluster objects for both 5' and 3' sides,
        keyed by cluster IDs.
    args : class
        Input arguments containing parameters such as the minimum anchor length
        and whether to include reference coordinates.
    temp_dir : str
        Path to the temporary directory where intermediate files are stored.

    Returns
    -------
    out : dict
        A dictionary containing the output for each anchor cluster, keyed by
        side ('5' or '3').

    """
    outhandle = open(os.path.join(args.outpath, f'{args.prefix}.anchors.tsv'), 'w')

    header = ['anchor_id', 'side', 'num_reads', 'consensus']
    
    if args.detailed:
        if args.reference:
            header += ['ref', 'ref_start', 'ref_end', 'ref_strand', 'ref_cigar', 'ref_mapq']
        header += ['target_start', 'target_end','anchor_slope', 'anchor_intercept', 'prop_sites_with_mismatches']
        
    outhandle.write('\t'.join(header) + '\n')

    for side in cluster_d:
 
        for cluster_id in cluster_d[side]:
            cl = cluster_d[side][cluster_id]
            target_pos = [i for i in cl.target_cov]
                       
            row = [cl.cluster_id, side, cl.nr_reads, str(cl.consensus.seq)]
            
            if args.detailed:
                if args.reference:
                    row += [cl.ref, cl.ref_start, cl.ref_end, cl.ref_strand, cl.ref_cigar, cl.ref_mapq]
                row += [min(target_pos), max(target_pos),  cl.anchor_lm[1], cl.anchor_lm[0], cl.prop_sites_with_mismatches]

            outhandle.write('\t'.join(map(str, row)) + '\n')
                        
    outhandle.close()


def write_reference_output(overlaps, cluster_d, args):
    
    """
    Write the reference output with insertion details to a file or stdout.

    This function processes identified insertions and writes details such as 
    chromosome, position, strand, support from anchor reads, and target site 
    duplication (TSD) to a specified output file or standard output. If provided, 
    it also includes gene information and the distance to the gene.

    Parameters
    ----------
    args : class
        Input arguments containing the output file path, annotation file, and reference file.
    overlaps : list
        List of identified insertions with details on position, cluster numbers, and TSD length.
    cluster_d : dict
        Dictionary containing clusters of anchor reads for both 5' and 3' sides.
    outpath : str
        Path to the output directory (not used in this function).

    Returns
    -------
    None
    """

    outhandle = open(os.path.join(args.outpath, f'{args.prefix}.reference_insertions.tsv'), 'w')
        
    header = ['chromosome', 'position', 'strand', 'TSD', 'support_5', 'support_3','anchor_5', 'anchor_3']
    
    # If an annotation is provided, load it and add gene information to output
    if args.annot:
        header += ['gene', 'dist_to_gene']
        
        annot = pandas.read_csv(
            args.annot, sep='\t', comment='#', 
            names=['seqid', 'source', 'type', 'start', 'end','score', 'strand', 'phase','attributes'])
        
        # Remove CDS entries
        annot = annot[annot['type'].isin(['gene', 'pseudogene', 'mobile_genetic_element'])]
        annot = annot.reset_index(drop=True)
        
    if args.detailed:
        header += ['mapq_5', 'mapq_3', 'cigar_5', 'cigar_3']
        
    # Get chromosome length
    reference = SeqIO.read(args.reference, 'fasta')
    chrom_length = len(reference.seq)
    outhandle.write('\t'.join(header) + '\n')
    
    # Get chromosome name, assuming that the reference is a single contig    
    chromosomes = [seq_record.id for seq_record in SeqIO.parse(args.reference, 'fasta')]
    chrom = chromosomes[0]

    # Now loop through identified mutations 
    for ins in overlaps:

        position = ins[1]
        strand = '+' if ins[0].startswith('5prime') else '-'

        # Nr anchor reads
        if strand == '+':
            five_cl_nr = int(ins[0].split('_')[1])
            three_cl_nr = int(ins[2].split('_')[1])
        elif strand == '-':
            five_cl_nr = int(ins[2].split('_')[1])
            three_cl_nr = int(ins[0].split('_')[1])
            
        support_5 = len(cluster_d['5'][five_cl_nr].reads)
        support_3 = len(cluster_d['3'][three_cl_nr].reads)

        # TSD
        tsd_len = ins[4]
        tsd = cluster_d['5'][five_cl_nr].consensus.seq[-tsd_len:]

        # Anchor ID, mapq and cigar
        anchor5_id = ins[0]
        anchor3_id = ins[2]
        mapq_5 = cluster_d['5'][five_cl_nr].ref_mapq
        cigar_5 = cluster_d['5'][five_cl_nr].ref_cigar
        mapq_3 = cluster_d['3'][three_cl_nr].ref_mapq
        cigar_3 = cluster_d['3'][three_cl_nr].ref_cigar

        outline = [chrom, str(position), strand, str(tsd), str(support_5), str(support_3),anchor5_id, anchor3_id]
        
        if args.annot:
            gene, dists_to_gene = gene_overlap(position, annot, chrom_length)
            outline += [gene, dists_to_gene]
            
        if args.detailed:
            outline += [str(mapq_5), str(mapq_3), cigar_5, cigar_3]

        outline = map(str, outline)

        outhandle.write('\t'.join(outline) + '\n')
        
    outhandle.close()
        
    
def remove_outliers(lista):
    """ Removes outliers from a list of integers, where outliers are defined
    as in boxplots
    
    Parameters
    ----------
    lista : list
        A list of integers.

    Returns
    -------
    lista_filt : list
        List of integers with outliers removed.

    """
    q_1 = numpy.percentile(lista, 25)
    q_3 = numpy.percentile(lista, 75)

    boxlength = q_3 - q_1

    lower = max(min(lista), q_1 - 1.5*boxlength)
    upper = min(max(lista), q_3 + 1.5*boxlength)

    lista_filt = [x for x in lista if x >= lower and x <= upper]
    return lista_filt



def consensus_from_bam(region, bamfile, min_mapq=1, min_baseq=1):
    """ Create a pileup file for a region and extract consensus sequence.
    Ignores indels.
    

    Parameters
    ----------
    region : list
        Genomic region with chromosome, start, end, 1-based.
    bamfile : str
        Path to bam file.
    min_mapq : int, optional
        Minimum mapping quality [1].
    min_baseq : int, optional
        Minimum base quality [1].
        
    Returns
    -------
    consensus_seq : str
        Consensus sequence of the region.

    """
    pybam = pysam.AlignmentFile(bamfile, "rb")

    chrmsm, strt, end = region[0], region[1]-1, region[2]

    pile_dict = dict()
    crap_reads = set()

    for pileupcolumn in pybam.pileup(chrmsm, strt, end, **{"truncate": True}):
        pos = pileupcolumn.pos
        pile_dict[pos] = []

        for pileupread in pileupcolumn.pileups:

            read = pileupread.alignment

            if read.query_name in crap_reads:
                continue

            if pileupread.is_del or pileupread.is_refskip:
                continue

            elif read.mapping_quality < min_mapq:
                crap_reads.add(read.query_name)
                continue

            elif read.query_qualities[pileupread.query_position] < min_baseq:
                continue

            base = pileupread.alignment.query_sequence[pileupread.query_position]
            pile_dict[pos].append(base)

    consensus_seq = ''
    for i in range(strt, end):
        try:
            bases = pile_dict[i]
            cov = len(bases)
            if cov == 0:
                consensus_seq += 'N'
                continue
            base_counts = Counter(bases)

            # most common base
            cons_base_list = base_counts.most_common(1)
            consensus_base = max(cons_base_list, key=lambda x: x[1])[0]

            consensus_seq += consensus_base

        except KeyError:
            consensus_seq += 'N'
            continue

    pybam.close()
    return consensus_seq
