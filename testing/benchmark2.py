#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Benchmark detettore using reads simulated from assembly

IS insertions are identified by creating a graph including the reference
and the assembly from which reads are simulated, and then mapping the IS
library against the graph
Created on Tue Apr 18 12:12:57 2023
@author: cristobal


GAF format (https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf)
PAF format (https://lh3.github.io/minimap2/minimap2.html)


class args:

    def __init__(self):

        self.REF = '/home/cristobal/TB/projects/detettoreTB/github/detettore6110/resources/reference/MTBC0_v1.1.fasta'
        self.ASM = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/A_ancestral_reference/genomes/NC_021251.1.fasta'
        self.IS_SEQ = '/home/cristobal/TB/projects/detettoreTB/github/detettore6110/resources/is_targets/IS6110.fasta'
        self.DETETTORE_RESULTS = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/C_detettore/simulation_results/CCDC5079/CCDC5079.Anc01_RL150_COV50.resultati.tsv'
        self.PREF = 'RW-TB008'

args = args()


"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import pandas

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_args():

    parser = argparse.ArgumentParser(
        description=
        'Compare the results of detettore6110 against the truth. This assumes \
            that reads were simulated from a complete assembly and used as input \
                for detettore. This script compares the reference against the \
                    assembly using minigraph and extract bubbles that correspond \
                        to IS polymorphisms.')

    parser.add_argument(
        '-ref', dest='REF',
        required=True,
        help='Reference genome')

    parser.add_argument(
        '-asm', dest='ASM',
        required=True,
        help='Assembled genome')

    parser.add_argument(
        '-is', dest='IS_SEQ',
        required=True,
        help='Target insertion sequence')

    parser.add_argument(
        '-res', dest='DETETTORE_RESULTS',
        required=True, nargs='+',
        help='Results of detettore6110.')

    args = parser.parse_args()

    return args


def align_assemblies(reference, assembly, min_len=0, sv_to_file=False):
    """ Align to genomes and call variants, using minimap2 and paftools.

    Args:
        reference (str): Path to reference genome
        assembly (str): Path to assembly from which reads were simulated
        min_len (int): Minimum length of variants
        svs_to_file (bool): Write structural variants to file. Default False.
    Returns:
        str: Path to the temporary fasta file containing SVs
        dict: dictionary with variants, with position as key and vcf line as value
        
    """

    variants_d = {}
    svs = []
    
    # Create a named temporary file
    fd, variants_fasta = tempfile.mkstemp(prefix='structural_variants', suffix='.fa')
    
    map_proc = subprocess.Popen(
        ('minimap2', '-c', '-x', 'asm5', '--cs', reference, assembly),
        stdout=subprocess.PIPE)
    
    sort_proc = subprocess.Popen(
        ('sort', '-k6,6', '-k8,8n'),
        stdin=map_proc.stdout, stdout=subprocess.PIPE)
    
    var_raw = subprocess.check_output(
        ('paftools.js', 'call', '-f', reference, '-'),
        stdin=sort_proc.stdout).splitlines()
    
    for line in var_raw:
        
        line = line.decode('utf-8')

        if line.startswith('#'):
            continue
    
        fields = line.strip().split('\t')
  
        ref = fields[3]
        alt = fields[4]
        pos = fields[1]
   
        if len(ref) == 1 and len(alt) == 1:
            continue

        var_len = abs(len(ref) - len(alt))
        if var_len < min_len:
            continue
        
        variants_d[pos] = fields
 
        typus = 'INS' if len(alt) > len(ref) else 'DEL'
 
        seq = ref if typus == 'DEL' else alt 
  
        var_id = fields[0] + '-' + fields[1] + '-' + typus

        rec = SeqRecord(
            Seq(seq), id=var_id, name=var_id, description='')
        
        svs.append(rec)
    
    with os.fdopen(fd, 'w') as temp_file:
        SeqIO.write(svs, temp_file, 'fasta')

    return variants_fasta, variants_d



def variants_vs_is(variants_fasta, is_fasta):
    """ Map IS sequence against structural variants. 
    
    Try blast instead of minimap??     
    
    
    Args:
        variants_fasta
        is_fasta
    Returns:
        
    """
    
    truth = {}
    map_proc = subprocess.Popen(
        ('minimap2', is_fasta, variants_fasta),
        stdout=subprocess.PIPE
        )
    
    map_out = map_proc.stdout.read().splitlines()
    
    for line in map_out:
                
        line = line.decode('utf-8')
        
        fields = line.strip().split()

        varname = fields[0].split('-')
        pos = int(varname[1])
 

  
        d = {
            'type': varname[2],
  
            'query_length': fields[1],
            'query_start': fields[2],
            'query_end': fields[3],
    
            'target_name': fields[5],
            'target_length': fields[6],
            'target_start': fields[7],
            'target_end': fields[8],

            'aln_block_length': fields[10],  # total nr of matches, mismatches, indels

            }

        truth[pos] = d

    return truth

        
def compare_results_with_truth(args, results, truth_d, imprec=100):
    """_summary_
    
    The "true" position is off a few bp because of the target site duplication:
        usually one of the duplicate sequence is included in the SV fasta 

    Args:
        results (_type_): _description_
        truth (_type_): _description_
        imprec (int): 
            
    """

    resultati = pandas.read_csv(args.DETETTORE_RESULTS, sep='\t')
    
    true_positives = {}
    others = {}
    
    performance = {
        'tp' : [],
        'fp' : [],
        'fn' : [],
        'doublecheck' : []
        }
    
    for pos in truth_d:
        
        if truth_d[pos]['type'] != 'INS':
            continue
         
        closest = (int, float('inf'))
        
        for i, row in resultati.iterrows():
            
            r_pos = row['position']
            
            tsd = row['TSD']
            tsd_l = 0 if pandas.isna(tsd) else len(row['TSD'])  
            pos_corr = r_pos - tsd_l
            dist = abs(pos - pos_corr)
            
            if dist < closest[1]:
                closest = (i, dist)

        print(closest)
                
        # True positive
        if closest[1] == 0:
            
            true_positives[pos] = closest
            
        else:
            others[pos] = closest
            

    # Indices of true positives in detettore results            
    tp_i = [true_positives[x][0] for x in true_positives]
            
    
    # False negatives: In truth but not in detettore results
    false_negatives = {}
    doublecheck = {} # 
    
    for pos in truth_d:
        
        if pos in true_positives:
            continue
        
        elif pos in others:
            
            i = others[pos][0]
            
            if i in tp_i:
                false_negatives[pos] = others[pos]
            
            else:
                doublecheck[pos] = others[pos]
                
    # False positives: in detettore results but not in truth
    false_positives = [i for i, row in resultati.iterrows() if i not in tp_i]
    
            
            
            
        
        
        



        

#%%

def main():

    args = get_args()

    # Compare assemblies and get variants
    variants_fasta, variants_d = align_assemblies(args.REF, args.ASM, min_len=1000)

    # Map variants against IS and get reference positions
    # of IS polymorphisms
    truth_d = variants_vs_is(variants_fasta, args.IS_SEQ)

    # Compare detettore positions with true reference positions
    performance = compare_results_with_truth(args, args.DETETTORE_RESULTS, truth_d)
    
    
    
    
    
    
    
    
    # Clean up
    os.remove(sv_fasta_path)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
