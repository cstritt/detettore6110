#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Module for detecting insertion sites based on mapping against a reference
"""


import re
import numpy

from collections import Counter


class ISpoly:
    """ IS polymorphism object in which all evidence is combined.
    """

    def __init__(self):

        self.chromosome = str()
        self.position = int()
    
        # Set with names of clustering reads
        self.reads = []
        self.te_hits = {}
        self.breakpoint = [ [] , [] ]
        self.region_start = float('inf')
        self.region_end = -float('inf')

        
    def evaluate_reads(self, reads, hits):
        """ 

        Parameters
        ----------
        reads : list
            A list of pysam reads clustering by position.
        hits : dict
            Dictionary with read IDs as keys and mapping of read against
            IS as value (output of parse_paf function).

        Returns
        -------
        Fills the ISpoly object above with all information required to 
        generate the output table.

        """
        
        for read in reads:
            
            cigar = read.cigartuples

            s = 1 if cigar[0][0] == 4 else 0

            # Region for local coverage
            if read.reference_start < self.region_start:
                self.region_start = read.reference_start
            if read.reference_end > self.region_end:
                self.region_end = read.reference_end
                
            # Insertion breakpoint
            break_pos = read.reference_end if s == 0 else read.reference_start
            self.breakpoint[s].append(break_pos)

            # TE hits dictionary
            hit = hits[read.query_name]
            target = hit['target_name']

            if target not in self.te_hits:

                self.te_hits[target] = {
                    'strand' : {'+' : 0, '-' : 0},
                    'hit_mapqs' : {},
                    'anchor_mapqs' : {},
                    'combined_mapqs' : {},
                    'aligned_positions_left' : set(),
                    'aligned_positions_right' : set(),
                    'side' : [0,0]
                    }

            # Counting support from down- and upstream
            self.te_hits[target]['side'][s] += 1

            strand = '+' if hit['strand'] == '+' else '-'

            self.te_hits[target]['strand'][strand] += 1

            self.te_hits[target]['hit_mapqs'][read.query_name] = hit['mapping_quality']
            self.te_hits[target]['anchor_mapqs'][read.query_name] = read.mapping_quality
            self.te_hits[target]['combined_mapqs'][read.query_name] = hit['mapping_quality'] + read.mapping_quality

            # Nr aligned positions on the target
            if s == 0:
                self.te_hits[target]['aligned_positions_left'].update(
                    range(hit['target_start'], hit['target_end'])
                    )
            elif s == 1:
                self.te_hits[target]['aligned_positions_right'].update(
                    range(hit['target_start'], hit['target_end'])
                    )


def parse_paf(paf_file, min_aln_len):
    """ Read paf output of minimap2 and put information into a dictionary
    
    If there are supplementary alignments, the mapping quality of a read is
    recalibrated. (Relevant if mapping qualities are used to calculated 
    genotype likelihoods.)
    

    Parameters
    ----------
    paf_file : str
        Path to paf file.
    min_aln_len : int
        Minimum alignment length.

    Returns
    -------
    minimapd : dict
        A dictionary with query names as keys and dictionaries with paf 
        information as values.

    """
    
    minimapd = {}
    
    with open(paf_file) as f:
        
        for line in f:
            
            fields = line.strip().split()
            
            if int(fields[10]) < min_aln_len:
                continue
            
            readname = fields[0]
        
            d = {
        
                'strand' : fields[4],
        
                'target_name' : fields[5],
                'target_start' : int(fields[7]),
                'target_end' : int(fields[8]),
                'target_range' : set(range(int(fields[7]), int(fields[8]))),
        
                'aln_block_length' : int(fields[10]), # total nr of matches, mismatches, indels
                'mapping_quality' : int(fields[11]),
        
                'tp' : re.search(r'tp:A:(.*?)\t', line).group(1),
                'cm' : int(re.search(r'cm:i:(.*?)\t', line).group(1)),
                's1' : int(re.search(r's1:i:(.*?)\t', line).group(1)),
                's2' : int(re.search(r's1:i:(.*?)\t', line).group(1))
        
                }
        
            # Secondary alignments: if they align to the same element, recalibrate mapping quality..
            if readname in minimapd:
        
                if d['target_name'] == minimapd[readname]['target_name']:
        
                    d_primary = minimapd[readname]
                    xmin = 1 if (minimapd[readname]['cm'] / 10) >= 1 else minimapd[readname]['cm']
                    mapq_recal = 40 * xmin * numpy.log(d_primary['s1'])
                    
                    if mapq_recal > 60:
                        mapq_recal = 60
                    
                    minimapd[readname]['mapping_quality'] = mapq_recal
        
                    # if s1 and s2 are equal, add positions form secondary alignment
                    if d['s1'] == d_primary['s1']:
                        minimapd[readname]['target_range'].update(d['target_range'])
        
            else:
                minimapd[readname] = d
    
    return minimapd


def TE_hit_summary(reads, hits):
    """ For a set of clustering reads, create a summary of their mapping 
    against the TE library.
    
    Relevant if target fasta contains multiple IS sequences...

    Parameters
    ----------
    reads : list
        A list of pysam reads clustering by position.
    hits : dict
        A dictoinary, output of parse_paf(), showing if and how a read
        maps against the IS target

    Returns
    -------
    te_hits : dict
        DESCRIPTION.
    best : str
        Name of the IS with the highest score (cumulative mapping quality)

    """
    
    te_hits = {}
    
    for read in reads:
        
        cigar = read.cigartuples

        s = 1 if cigar[0][0] == 4 else 0

        # TE hits dictionary
        hit = hits[read.query_name]
        target = hit['target_name']

        if target not in te_hits:

            te_hits[target] = {
                'strand' : {'+' : 0, '-' : 0},
                'hit_mapqs' : {},
                'anchor_mapqs' : {},
                'combined_mapqs' : {},
                'aligned_positions' : set(),
                'side' : [0,0]
                }

        # Counting support from down- and upstream
        te_hits[target]['side'][s] += 1


        strand = '+' if hit['strand'] == '+' else '-'

        te_hits[target]['strand'][strand] += 1

        te_hits[target]['hit_mapqs'][read.query_name] = hit['mapping_quality']
        te_hits[target]['anchor_mapqs'][read.query_name] = read.mapping_quality
        te_hits[target]['combined_mapqs'][read.query_name] = hit['mapping_quality'] + read.mapping_quality

        # Nr aligned positions on the target
        te_hits[target]['aligned_positions'].update(
            range(hit['target_start'], hit['target_end'])
            )
        
    # Get highest scorting TE
    score = -1
    best = ''

    for te in te_hits:
        mapq_score = sum([te_hits[te]['hit_mapqs'][x] for x in te_hits[te]['hit_mapqs']])
        if mapq_score > score:
            score = mapq_score
            best = te
            
    return te_hits, best
    

def cluster_splitreads(splitreads, hits, MIN_SPLIT_COV=5):
    """ Cluster splitreads based on the exact positions at which reads
    are cut.
    

    Parameters
    ----------
    splitreads : dict
        DESCRIPTION.
    hits : dict
        A dictoinary, output of parse_paf(), showing if and how a read
        maps against the IS target.
    MIN_SPLIT_COV : int, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    split_positions : dict
        A dictionary with (position, side) as keys, a list of clustering
        reads as values.

    """
    
    split_positions = {}

    # Discard splitreads not mapping to IS
    split_spurious = set(splitreads.keys()) -  set(hits.keys())

    for readname in split_spurious: 
        del splitreads[readname]
        
        
    # dictionary of split positions
    for readname in splitreads:
        
        for read in splitreads[readname]:
        
            cigar = read.cigartuples
            
            if len(cigar) > 1: 
            
                first, last = cigar[0], cigar[-1]
                
                if first[0] in (4,5):
                    pos = read.reference_start

                    side = "right"
      
                elif last[0] in(4,5):
                    pos = read.reference_end             
                    side = "left"
                    
                k = (pos, side)
                    
                if k not in split_positions:
                    split_positions[k] = []
                        
                split_positions[k].append(read)
                
                
    remove = []

    for k in split_positions:
        if len(split_positions[k]) < MIN_SPLIT_COV:
                remove.append(k)
                
    for k in remove:
        split_positions.pop(k, None)           
    
    return split_positions

    
def merge_clusters(split_positions, hits, tsd_max_len = 10, max_dist = 10000):
    """ 
    Go through clusters of split reads and merge them when they point to 
    the same insertion event.
    
    In case of precise insertions, L and R clusters will overlap. Imprecise 
    insertion sites arise when SVs have messed up the insertion signature.
    Often this results in L and R clusters that are separated, e.g.:
    
    MTBC0	1993886	-		L	16	0	16	1258	1353
    MTBC0	1994193	-		R	0	8	8	1	81
    
    
    Parameters
    ----------
    split_positions : dict
        DESCRIPTION.
    hits : dict
        A dictoinary, output of parse_paf(), showing if and how a read
        maps against the IS target. Only reads are kept which are present
        in this dictionary.        
    tsd_max_len : int, optional
        Maximum lenght of the target site duplication. The default is 10.
    max_dist : int, optional
        Maximum distance between L and R support for an insertion event. 
        The default is 10000.

    Returns
    -------
    clusters : list
        A list with insertion event condidates, where each element is a tuple 
        (position, evidence type, ISpoly). Evidence type is L, R, LR, or 
        LR_imprecise. 
        
        At present, detettore6110 is designed to output only clear LR 
        signatures (maximize specificity).

    """
    
    pos_sorted = sorted(split_positions.keys())

    clusters = []
    skip = False 

    for i in range(0, len(pos_sorted) - 1):
        
        if skip:
            skip = False
            continue
        
        pos = pos_sorted[i]
        reads = split_positions[pos]
        te_hits, best = TE_hit_summary(reads, hits)
        
        next_pos = pos_sorted[i+1]
        dist = next_pos[0] - pos[0]
        
        
        # Precise insertion with TSD
        if pos[1] == "right" and next_pos[1] == "left" and dist <= tsd_max_len:
            reads = split_positions[pos] + split_positions[next_pos]
            cat = "LR"
            skip = True
            
        # Imprecise: L and R clusters with large or no overlap
        # Assume same insertion event if same strand and opposite ends of the IS
        # are covered
        
        elif (pos[1] == "right" and next_pos[1] == "left" and dist > tsd_max_len) or \
            (pos[1] == "left" and next_pos[1] == "right" and dist < max_dist):
                
            if is_same_event(
                    (split_positions[pos], split_positions[next_pos]), hits
                    ):
                
                reads = split_positions[pos] + split_positions[next_pos]
                cat = "LR_imprecise"
                skip = True
                
            else:
                reads = split_positions[pos]
                cat = "L" if pos[1] == "left" else "R"
                
        # Only one side (because the other is in a sequence not present in the reference)
        else:
            reads = split_positions[pos]
            cat = "L" if pos[1] == "left" else "R"
        
        cl = ISpoly()
        cl.evaluate_reads(reads, hits)
        
        clusters.append((pos[0],cat, cl))
        
    return clusters


def is_same_event(read_clusters, hits):
    """
    Check if two clusters of reads support the same insertion event. 
    
    Criteria:
        - both clusters support the same strand
        - reads of the two clusters reach into opposite ends of the IS
    
    
    Parameters
    ----------
    reads : tuple
        Tuple of length two with read clusters from both sides.
    hits : dict
        Output of parse_paf(), showing if and how a read
        maps against the IS target. Only reads are kept which are present
        in this dictionary..

    Returns
    -------
    True if the L and the R cluster have reads that reach into the opposite
    ends of the IS and support the same insertion strand.

    """
    
    strand = {
        0 : [],
        1 : []
        }
    
    aligned_positions = {
        0 : set(),
        1 : set()
        }
    
    for i, cluster in enumerate(read_clusters):
        
        for read in cluster:
            
            hit = hits[read.query_name]
           
            s = '+' if hit['strand'] == '+' else '-'
            strand[i].append(s)
  
            
            aligned_positions[i].update(
                range(hit['target_start'], hit['target_end'])
                )
            
    strand0 = Counter(strand[0]).most_common(1)[0][0]
    strand1 = Counter(strand[1]).most_common(1)[0][0]

    aln_overlap = aligned_positions[0].intersection(aligned_positions[1])
    
    if strand0 != strand1 or aln_overlap:
        return False
    else:
        return True
    