#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Module for writing output

"""

import sys
import numpy
import pysam

from Bio import SeqIO
from collections import Counter


def write_output(args, clusters, copy_number, outpath, 
                 both_sides=True, require_tsd=True, mapq_filt=False):
    """ Write detettore6110 output, including copy number and 
    insertion sites. 
    
    Bug when required_tsd=False:
        still tries to get tsd, running into position error when start > end
           
        
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
        
        'TSD', 'position_L', 'position_R', 
    
        'IS_start', 'IS_end'
        
        ]
    
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
        
        TE_L_START = min(te_hits['aligned_positions_left'])
        TE_L_END = max(te_hits['aligned_positions_left'])
        TE_START = f'{TE_L_START},{TE_L_END}'
        
        TE_R_START = min(te_hits['aligned_positions_right'])
        TE_R_END = max(te_hits['aligned_positions_right'])
        TE_END = f'{TE_R_START}, {TE_R_END}'
        
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
            TSD, POS_LEFT, POS_RIGHT, 
            TE_START, TE_END
                   ]        
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
