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
    
    
def exit_handler(args, temp_dir, working_dir):
    """ Cleanup after program finish. If --keep is given, copy contents of 
    temporary directory to working directory before deleting it"""
    if args.keep:  # Copy contents of temporary to working directory
        shutil.copytree(temp_dir, os.path.join(working_dir, os.path.basename(temp_dir)))
        
    shutil.rmtree(temp_dir)
    


def gene_overlap(position, annotation):
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

    Returns
    -------
    A list with one element for each position, containing the gene name 
    and the distance to the gene. If the position is in a gene, this is one
    name and distance 0; if the position is between genes, this is two genes 
    and two distances separated by ;   

    """
    
    regex_patterns = [
        r";gene=([^;\n]+)", 
        r"locus_tag=([^;\n]+)", 
        r"mobile_element_type=([^;\n]+)"
        ]
    
    # Find the closest intervals on either side of each position
    idx_s = bisect.bisect_right(annotation['start'], position)
    idx_e = bisect.bisect_right(annotation['end'], position)
    
    # Overlapping feature
    if idx_e == idx_s - 1:
                
        gene_info = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_e])
        
        return [gene_info, '0']
                                
    # Inbetween features
    elif idx_s == idx_e:
        
        gene_info_5 = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_s - 1])
        
        dist_to_5 = position - annotation['end'][idx_s - 1]
        
        gene_info_3 = gene_id_regex(
            regex_patterns, annotation['attributes'][idx_s])
        
        dist_to_3 = annotation['start'][idx_s] - position
        
        return [f'{gene_info_5};{gene_info_3}', f'{dist_to_5};{dist_to_3}']
        

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
        

    

def write_output(args, clusters, copy_number, outpath, 
                 both_sides=True, require_tsd=True, mapq_filt=True):
    """ Write detettore6110 output, including copy number and 
    insertion sites. 
    
    (Only the mapq filter works properly at the moment! All other need to
    be set to true. Maybe implement less stringent filtering and allow
    less clear-cut insertion signatures, with separated or single cluster
    support.)
           
        
    Parameters
    ----------
    args : args
        Arguments passed to detettore.
    clusters : list
        Split read clusters obtained from cluster_splitreads function.
    copy_number : int
        IS copy number from copynumbers.get_copy_number.
    outpath : str
        Write file to this directory.
        
    both_sides : bool, optional
        Require split reads from 5' and 3' side. Default true. 
    require_tsd : bool, optional
        Require that splitread clusters overlap because of target site 
        duplication. Default true.
    mapq_filt : bool, optional
        Require that anchor read mapq do not sum to 0. Default false.
        
    
    Returns
    -------
    Writes detettore6110 results to stdout or to file.

    """

    if args.outfile:
        outhandle = open(args.outfile, 'w')
        
    header = [
        
        'chromosome', 'position', 'strand', 
        
        'support_L', 'support_R', 'mean_mapq',
        
        'TSD', 'position_L', 'position_R' 
    
        #'IS_start', 'IS_end'
        
        ]
    
    # If an annotation is provided, load it and add gene information to output
    if args.annot:
        header += ['gene', 'dist_to_gene']
        
        annot = pandas.read_csv(
            args.annot, sep='\t', comment='#', 
            names=['seqid', 'source', 'type', 'start', 'end','score', 'strand', 'phase','attributes'])
        
        # Remove CDS entries
        annot = annot[annot['type'].isin(['gene', 'pseudogene', 'mobile_genetic_element'])]
        annot = annot.reset_index(drop=True)
    
    headerstr = '\t'.join(header)
    firstlines = f'#CN {copy_number}\n{headerstr}\n'
    outhandle.write(firstlines) if args.outfile else sys.stdout.write(firstlines)
    

    # Get chromosome name, assuming that the reference is a single contig    
    chromosomes = [seq_record.id for seq_record in SeqIO.parse(args.ref, 'fasta')]
    CHROM = chromosomes[0]

    for i, c in enumerate(clusters):
        
        # Nr supporting reads
        SUPPORT_LEFT = len(c[2].breakpoint[0])
        SUPPORT_RIGHT = len(c[2].breakpoint[1])
        
        if both_sides and SUPPORT_LEFT == 0 or SUPPORT_RIGHT == 0:
                continue
        
            
        # IS info
        is_name = list(c[2].te_hits.keys())[0]
        te_hits = c[2].te_hits[is_name]
        
        #TE_L_START = min(te_hits['aligned_positions_left'])
        #TE_L_END = max(te_hits['aligned_positions_left'])
        #TE_START = f'{TE_L_START}-{TE_L_END}'
        
        #TE_R_START = min(te_hits['aligned_positions_right'])
        #TE_R_END = max(te_hits['aligned_positions_right'])
        #TE_END = f'{TE_R_START}-{TE_R_END}'
        
        STRAND = '+' if te_hits['strand']['+'] > te_hits['strand']['-'] else '-'
        
        
        # Position: the right-most base of the left cluster (5' strand) is considered the insertion site
        POS_LEFT = max(remove_outliers(c[2].breakpoint[0])) 
        POS_RIGHT = min(remove_outliers(c[2].breakpoint[1])) 
        
        LR_DIST = POS_RIGHT - POS_LEFT
        
        if require_tsd and (LR_DIST > 0) or (LR_DIST < -args.max_tsd_len):
                continue

        POS = POS_LEFT  
        
        
        # Mean mapq per read
        sum_anchor_mapqs = sum([te_hits['anchor_mapqs'][k] for k in te_hits['anchor_mapqs']])
        MEAN_MAPQ = int(sum_anchor_mapqs / (SUPPORT_LEFT + SUPPORT_RIGHT))


        # At least some anchors should have a mapq > 0
        if mapq_filt and sum_anchor_mapqs == 0:
                continue
        
        
        # Target site duplication
        region = [CHROM, POS_RIGHT + 1, POS_LEFT]
        TSD = consensus_from_bam(region, f'{outpath}/reads_vs_ref.bam')
        
        
        # Write output
        outline = [
            CHROM, POS, STRAND, 
            SUPPORT_LEFT, SUPPORT_RIGHT, MEAN_MAPQ,
            TSD, POS_LEFT, POS_RIGHT 
            #TE_START, TE_END
            ]
        
        # Add gene information
        if args.annot:
            gene_info = gene_overlap(POS, annot)
            try:
                outline += gene_info
            except TypeError:
                print(POS)

        outline = map(str, outline)

        if args.outfile:
            outhandle.write('\t'.join(outline) + '\n')
        
        else:
            sys.stdout.write('\t'.join(outline) + '\n')
    
    if args.outfile:    
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
