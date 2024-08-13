#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Benchmarking using a second high-quality assembly

Created on Tue Jun 15 15:50:30 2021
@author: cristobal


PAF format (https://lh3.github.io/minimap2/minimap2.html)



# Test run TB ASM1488462v1
class args:
    def __init__(self):
        
        self.reference = '/home/cristobal/TB/analyses/detettore/INPUT/MTB_ancestor_reference.fasta'
        self.assembly = '/home/cristobal/TB/analyses/cactus/genomes/PRJNA555636_Modlin/ASM1488462v1.fna'
        self.detettore_vcf = '/home/cristobal/TB/analyses/detettore/benchmark/simReads_Modlin/detettore/ASM1488462v1.detettore.vcf.gz'
        self.TElib = '/home/cristobal/TB/analyses/detettore/INPUT/IS_library.mycobacteria80perc.fa'
        self.keep = True
        
args = args()

# Test run Bdis
class args:
    def __init__(self):
        
        self.reference = '/home/cristobal/Bdis/data/reference_genome/v3.1/Bdistachyon_314_v3.0.fa'
        self.assembly = '/home/cristobal/Bdis/data/BdTR7a_assembly/Bdistachyon_BdTR7a.v1.0.fa'
        self.detettore_vcf = '/home/cristobal/github/detettore/benchmark/bdis/BdTR7a.detettore.vcf.gz'
        self.TElib = '/home/cristobal/databases/repeats/TREP/TREP_consensus_autonomous_simplifiedNames.fasta'
        self.keep = True
        
args = args()


"""

import argparse
import os
import re
import gzip
import pandas
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_args():

    parser = argparse.ArgumentParser(
        description='Detect TE polymorphisms from paired-end read data')

    parser.add_argument(
        "-r",
        dest="reference",
        required=True,
        help='Reference genome in fasta format.')

    parser.add_argument(
        '-a', 
        dest="assembly",
        required=True,
        help='Second high quality assembly.')
    
    parser.add_argument(
        '-t', 
        dest="TElib",
        required=True,
        help='TE consensus library, same as used with detettore.')
    
    parser.add_argument(
        "-v",
        dest="detettore_vcf",
        required=True,
        help='TIPs and TAPs called with detettore')
    
    parser.add_argument(
        '--keep',
        action='store_true',
        help='Keep intermediate files stored in tmp folder.')

    args=parser.parse_args()

    return args



def blastn(query, target):

    cmd = [
        'blastn',
        '-query', query,
        '-subject', target,
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send pident length',
        '-max_target_seqs', '1']

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = proc.stdout.read()
    return output.splitlines()




#%% MAIN
def main():
    
    args = get_args()  
    os.mkdir('tmp')
    
    #%% Load detettore results
    
    #detettore = pandas.read_table(args.detettore_vcf, comment='#' )
    
    detettore = {}
    
    with gzip.open(args.detettore_vcf) as f:
        for line in f:
            line = line.decode('utf-8')
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            # # skip invariant
            if '0/0' in fields[9]:
                continue
            
            if fields[0] not in detettore:
                detettore[fields[0]] = {}
                
            detettore[fields[0]][int(fields[1])] = fields




    #%% Test args
    
    
    class args:

        def __init__(self):

            self.REF = "example/MTB_ancestor_reference.fasta"
            self.ASM = "benchmark/genomes/N1202.fasta"
            self.TELIB = "example/isescan.is_seqs.fna.consensus_seqs.renamed.fa"
           

    args = args()  

        
    #%% Minimap: Align assemblies and call variants
    
    cmd = 'minimap2 -cx asm5 --cs %s %s | sort -k6,6 -k8,8n \
        | paftools.js call -f %s - > tmp/the_truth.vcf' % (args.REF, args.ASM, args.REF) 
        
    os.system(cmd)
    
    variants = []
    
    min_len = 50
    
    with open('tmp/the_truth.vcf') as f:
        for line in f:
            if line.startswith('#'):
                continue
    
            fields = line.strip().split('\t')
            
            ref = fields[3]
            alt = fields[4]
            
            if len(ref) == 1 and len(alt) == 1:
                continue
            
            var_len = abs(len(ref) - len(alt))
            if var_len < min_len:
                continue
            
            typus = 'INS' if len(alt) > len(ref) else 'DEL'
        
            seq = ref if typus == 'DEL' else alt 
            
            var_id = fields[0] + '-' + fields[1] + '-' + typus
            
            rec = SeqRecord(
                
                Seq(seq), id = var_id, name = var_id, description=''
                
                )
            
            variants.append(rec)
            
            
    SeqIO.write(variants, 'tmp/the_truth.min%ibp.fa' % (min_len) , 'fasta')
    
    
    #%% Map structural variants against TE library
    
    #blastout = blastn('tmp/the_truth.min%ibp.fa' % (min_len), args.TElib)
    
    cmd = 'minimap2 -xsr %s tmp/the_truth.min%ibp.fa \
        > tmp/the_truth-VS-TElib.paf' % (args.TELIB, min_len)
        
    os.system(cmd)
    
    truth = {}
    
    with open('tmp/the_truth-VS-TElib.paf') as f:
        for line in f:
            
            fields = line.strip().split()
        
            varname = fields[0].split('-')
            chromosome, pos = varname[0], int(varname[1])
            
            if chromosome not in truth:
                truth[chromosome] = {}
                
            d = {
                'type' : varname[2],
                
                'query_length' : fields[1],
                'query_start' : fields[2],
                'query_end' : fields[3],
                      
                'target_name' : fields[5],
                'target_length' : fields[6],
                'target_start' : fields[7],
                'target_end' : fields[8],
            
                'aln_block_length' : fields[10], # total nr of matches, mismatches, indels
       
                }
            
            truth[chromosome][pos] = d
        
        
    #%% False and true positives
    
    """ 
    Create 4 tables: TIPs 
    
    one with true and false positives,
    one with false negatives. The first contains all the detettore information
    and is used for Gaussian mixture modelling, the second contains the 
    WGA information
    
    """
    
    tip_header = [
        'chr', 'pos', 'type', 'meinfo', 'me', 'DPADJ', 
        'GT', 'GQ', 'ref_supp', 'alt_supp', 
        
        'imprec', 'DR', 'SR', 'AL', 'AP', 'BPIQR',

        'wga_pos', 'var_len', 'var_target', 'var_target_len', 'aln_len', 
        'category'
        ]
    
    tap_header = [
        'chr', 'pos', 'type', 'meinfo', 'me', 'DPADJ', 
        'GT', 'GQ', 'ref_supp', 'alt_supp', 
        
        'PDEV',

        'wga_pos', 'var_len', 'var_target', 'var_target_len', 'aln_len', 
        'category'
        ]
    
    positives = {
        'INS' : [],
        'DEL' :[]
        }
    
    # Allow some positional uncertainty
    pos_range = 20
    
    # Store true positives to identifiy false negatives
    tp = []
    
    
    for chromosome in detettore:
        
        try:
            true_positions = list(truth[chromosome].keys())
        except KeyError:
            print("Chromosome not in truth set:" + chromosome + '\n')
            continue
        
        for pos in detettore[chromosome]:
            
            fields = detettore[chromosome][pos]
           
            # Parse VCF line
            chromosome, pos = fields[0], int(fields[1])
            meinfo = re.search(r'MEINFO=(.*?);', fields[7]).group(1)
            me = meinfo.split(',')[0]
            svtype = re.search(r'SVTYPE=(.*?);', fields[7]).group(1)
            dpadj = re.search(r'DPADJ=(.*?);', fields[7]).group(1)
            
            outline = [chromosome, pos, meinfo, me, svtype, int(dpadj)]
            
            
            # GT fields 
            gt_pre = fields[9].split(':')
            gt = gt_pre[0]        
            gq = gt_pre[1]
            ref_supp = gt_pre[2].split(',')[0]
            alt_supp = gt_pre[2].split(',')[1]
            
            outline += [gt, int(gq), int(ref_supp), int(alt_supp)]
            

            # INFO fields
            if svtype == 'INS':
                imprec = '1' if 'IMPRECISE' in fields[7] else '0'
                dr = re.search(r'DR=(.*?);', fields[7]).group(1)
                sr = re.search(r'SR=(.*?);', fields[7]).group(1)
                al = re.search(r'AL=(.*?);', fields[7]).group(1)
                ap = re.search(r'AP=(.*?);', fields[7]).group(1)
                bpiqr = re.search(r'BPIQR=(.*?)$', fields[7]).group(1)
                
                outline += [imprec, int(dr), int(sr), int(al), float(ap), int(bpiqr)]
                
                
            if svtype == 'DEL':
                pdev = re.search(r'PDEV=(.*?)$', fields[7]).group(1)
                outline += [float(pdev)]
                
    

            # Check if present in truth
            closest = min(true_positions, key=lambda x:abs(x-pos))
            
            if abs(closest-pos) < pos_range:
                
                d = truth[chromosome][closest]
                wga_pos = closest
                var_len = int(d['query_length'])
                var_target = d['target_name']
                var_target_len = int(d['target_length'])
                aln_len = int(d['aln_block_length'])
                cat = 'TP'
                
                tp.append((chromosome, closest))
                  
            else:
                wga_pos = closest
                var_len = 'NA'
                var_target = 'NA'
                var_target_len = 'NA'
                aln_len = 'NA'
                cat = 'FP'
                
            outline += [str(wga_pos), var_len, var_target, var_target_len, aln_len, cat]
            positives[svtype].append(outline)
  
            
  
    
positives['INS'] = pandas.DataFrame(positives['INS'], columns = tip_header)
positives['DEL'] = pandas.DataFrame(positives['DEL'], columns= tap_header)



     
    #%% False negatives 
    n_header = ['chromosome', 'pos', 'var_len', 'te_name', 'te_len', 'aln_len']
    
    negatives = [ n_header ]
    
    for chromosome in truth:
        for pos in truth[chromosome]:
            
            if (chromosome, pos) not in tp:
                
                d = truth[chromosome][pos]
                 
                outline = [
                    chromosome,
                    str(pos),
                    d['query_length'],
                    d['target_name'],
                    d['target_length'],
                    d['aln_block_length']
                    ]
                
                negatives.append(outline)
            
            
    #%% Write tables
    with open('benchmark_positives.tsv', 'w') as f:
        f.write('\t'.join(p_header) + '\n')
        for loc in positives:
            f.write('\t'.join(positives[loc]) + '\n')
            
    with open('benchmark_negatives.tsv', 'w') as f:
    
        for loc in negatives:
            f.write('\t'.join(loc) + '\n')
                
    if not args.keep:
        os.rmdir('tmp')
        
    
if __name__ == '__main__':
    main()
      
        

# #%% Minigraph: create graph and call variants

# minigraph = 'minigraph -xggs -t4 %s %s > minigraph/mg.detettore_bm.gfa' % (ref_p, asm_p)
  
# gfa2fasta = 'gfatools gfa2fa -s minigraph/mg.detettore_bm.gfa \
#     > minigraph/mg.detettore_bm.fa'
  
# callvars = 'gfatools bubble minigraph/mg.detettore_bm.gfa > minigraph/mg.detettore_bm.var.bed'

# getASMpath = 'minigraph -xasm --call minigraph/mg.detettore_bm.gfa %s \
#     > minigraph/mg.detettore_bm.var.ASM.bed' % (asm_p)
        
        
# for cmd in [minigraph, callvars, getASMpath]:
#     os.system(cmd)
