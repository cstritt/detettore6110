#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Benchmark detettore using reads simulated from assembly

IS insertions are identified by creating a graph including the reference
and the assembly from which reads are simulated, and then mapping the IS
library against the graph.   
    

Created on Tue Apr 18 12:12:57 2023
@author: cristobal


GAF format (https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf)

1 	string 	Query sequence name
2 	int 	Query sequence length
3 	int 	Query start (0-based; closed)
4 	int 	Query end (0-based; open)
5 	char 	Strand relative to the path: "+" or "-"
6 	string 	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
7 	int 	Path length
8 	int 	Start position on the path (0-based)
9 	int 	End position on the path (0-based)
10 	int 	Number of residue matches
11 	int 	Alignment block length
12 	int 	Mapping quality (0-255; 255 for missing)


class args:

    def __init__(self):

        self.REF = '/home/cristobal/TB/projects/detettoreTB/github/detettore6110/resources/reference/MTBC0_v1.1.fasta'
        self.ASM = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/A_ancestral_reference/genomes/NC_021251.1.fasta'
        self.IS_SEQ = '/home/cristobal/TB/projects/detettoreTB/github/detettore6110/resources/is_targets/IS6110.fasta'
        self.PREF = 'RW-TB008'
       
args = args()


"""

import argparse
import sys
import subprocess
import re
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




#%% 

def get_args():

    parser = argparse.ArgumentParser(
        description='Compare the results of detettore6110 against the truth. \
            This assumes that reads were simulated from a complete assembly and used as input for detettore. \
                This script compares the reference against the assembly using minigraph and \
                    extract bubbles that correspond to IS polymorphisms.')
    
    parser.add_argument(
        '-ref', dest = 'REF',
        required = True,
        help = 'Reference genome')
    
    parser.add_argument(
        '-asm', dest = 'ASM',
        required = True,
        help = 'Assembled genome')
    
    parser.add_argument(
        '-is', dest = 'IS_SEQ',
        required = True,
        help = 'Target insertion sequence')
    
    parser.add_argument(
        '-res', dest = 'DETETTORE_RESULTS',
        required = True, nargs = '+',
        help = 'Results of detettore6110.')
    
    
    
    args=parser.parse_args()
    return args
    

def create_graph(ref, asm, telib, threads=4):
    """
    4) # GFA segments in the bubble including the source and the sink of the bubble
    (5) # all possible paths through the bubble (not all paths present in input samples)
    (6) 1 if the bubble involves an inversion; 0 otherwise
    (7) length of the shortest path (i.e. allele) through the bubble
    (8) length of the longest path/allele through the bubble
    (9-11) please ignore
    (12) list of segments in the bubble; first for the source and last for the sink
    (13) sequence of the shortest path (* if zero length)
    (14) sequence of the longest path (NB: it may not be present in the input samples)
    
    """
    create = ['minigraph', '-c', '-x','ggs', '-t', str(threads), ref, asm]

    with open('graph.gfa', "w") as outfile:
        subprocess.run(create, stdout=outfile, stderr=subprocess.DEVNULL)
        
        
    bubble_cmd = ['gfatools', 'bubble', 'graph.gfa']
    bubble_proc = subprocess.Popen(bubble_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    bubble_raw = bubble_proc.stdout.read().splitlines()
    variants = []
    for line in bubble_raw:
        line = line.decode('utf-8')
        fields = line.strip().split('\t')
        variants.append(fields)

    return variants 


#%% MAIN

def main():
    
    args = get_args()

    variants = create_graph(args.REF, args.ASM, args.IS_SEQ)
    
    
    # Write variant seqs to fasta
    seqs = []
    
    for v in variants:
        
        seq = Seq(v[13])
        name = '_'.join(v[:3])
        rec = SeqRecord(seq, id=name, name=name, description='')
        seqs.append(rec)
        
    SeqIO.write(seqs, 'bubbles.fasta', 'fasta')
    
    
    #%% Map against IS library
    
    cmd = ['minimap2', args.IS_SEQ, 'bubbles.fasta']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output_raw = proc.stdout.read().splitlines()
    
    minimapd = {}
    
    for line in output_raw:
    
        line = line.decode('utf-8')     
        fields = line.strip().split()
        
        var_id = fields[0]
    
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
        
        minimapd[var_id] = d
        
    for v in variants:
        name = '_'.join(v[:3])
        if name in minimapd:
            minimapd[name]['path'] = v[11]
            minimapd[name]['len_longest_path'] = v[7]
            minimapd[name]['len_shortest_path'] = v[6]
    
    
    #%% Output table
    
    header = ["chrom", "start", "end", "type", "IS", "path", "len_shortest_path", "len_longest_path", "aln_len"]
  
    sys.stdout.write('\t'.join(header) + '\n')
    
    for v in minimapd:
        
        region = v.split('_')
        chrom = region[0]
        start = region[1]
        end = region[2]
        
        typus = "TIP" if start == end else "TAP"
        
        outline = [
            chrom,
            start,
            end,
            typus,
            minimapd[v]['target_name'],
            minimapd[v]['path'],
            minimapd[v]['len_shortest_path'],
            minimapd[v]['len_longest_path'],
            str(minimapd[v]['aln_block_length'])
            ]
                
        sys.stdout.write('\t'.join(outline) + '\n')
        
    #%% Compare to truth
    
    
    
    
    
    
    
    
    
    
    
    
        
    #%% Clean up
    #os.remove('graph.gfa')
    os.remove('bubbles.fasta')
    
#%%
if __name__ == '__main__':
    main()
