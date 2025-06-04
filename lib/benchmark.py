# !/usr/bin/env python3

""" Benchmark detettore6110 with simulated reads

Given a set of assemblies annotated with ISEScan, extract the insertions matching
the target IS and create detettore like output in order to compare true insertions with 
those inferred from short reads (false positives, false negatives).

Add an assembly validation step: look for inverted repeats of the target IS outside annotated IS.

"""
#%%

import argparse
import os
import pandas
import sys
import subprocess

from Bio import Align
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_args():
    
    parser = argparse.ArgumentParser(
        
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    
        description="""
        Evaluate detettore6110 performance by comparing results from simulated reads with ISEScan annotation.
        Outputs true copy number, true positives, false negatives, false positives for both the 5' and the 3'
        anchor sequences
        """
    )

    #parser.add_argument('--is_target', help='Path to the target IS sequence in FASTA format')
    #parser.add_argument('--is_target_cluster', help='The cluster (subfamily) to which the target IS belongs according to ISEScan')
    parser.add_argument('--isescan', help='Path to ISEScan fasta output (.is.fna)')
    parser.add_argument('--assembly', help='Path to the assembly from which reads were simulated')
    parser.add_argument('--detettore', help='Path to the detettore anchors file (.anchors.tsv)')
    parser.add_argument('--min_id', type=float, default=0.9, help='Minimum identity threshold for matching IS and anchors')
    parser.add_argument('--anchor_length', type=int, default=300, help='Length of the true anchors extracted from the assembly')
    
    args = parser.parse_args()
    return args


def sequence_similarity(target, query):
    """
    Compute the similarity between a sequence and the target IS sequence.

    Parameters
    ----------
    seq : Seq or str
        The sequence to compare with the target IS sequence.
    is_target : Seq or str
        The target IS sequence.

    Returns
    -------
    similarity : float
        The similarity between the two sequences, computed as 1 - (number of mismatches / length of the alignment).

    """
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(target, query)
    # Get the best alignment
    best_alignment = alignments[0]
    # Count the number of mismatches
    mismatches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a != b)
    # Compute the distance as the number of mismatches divided by the length of the alignment
    distance = mismatches / len(best_alignment.target)
    similarity = 1 - distance
    return round(similarity, 2)


class true_insertions:
    def __init__(self):
        """
        Initialize a true_insertions object with target IS sequence, cluster, and minimum identity threshold.

        Parameters
        ----------
        is_target : str
            The path to the FASTA file containing the target IS sequence.
        is_target_cluster : str
            The cluster identifier for the target IS sequence.
        min_id : float
            The minimum identity threshold to consider a sequence as a true insertion.
        """
        
        target_seq = """
        TGAACCGCCCCGGCATGTCCGGAGACTCCAGTTCTTGGAAAGGATGGGGTCATGTCAGGTGGTTCATCGAGGAGGTACCCGCCGGAGCTGCGTGAGCGGGCG
        GTGCGGATGGTCGCAGAGATCCGCGGTCAGCACGATTCGGAGTGGGCAGCGATCAGTGAGGTCGCCCGTCTACTTGGTGTTGGCTGCGCGGAGACGGTGCGT
        AAGTGGGTGCGCCAGGCGCAGGTCGATGCCGGCGCACGGCCCGGGACCACGACCGAAGAATCCGCTGAGCTGAAGCGCTTGCGGCGGGACAACGCCGAATTG
        CGAAGGGCGAACGCGATTTTAAAGACCGCGTCGGCTTTCTTCGCGGCCGAGCTCGACCGGCCAGCACGCTAATTACCCGGTTCATCGCCGATCATCAGGGCC
        ACCGCGAGGGCCCCGATGGTTTGCGGTGGGGTGTCGAGTCGATCTGCACACAGCTGACCGAGCTGGGTGTGCCGATCGCCCCATCGACCTACTACGACCACA
        TCAACCGGGAGCCCAGCCGCCGCGAGCTGCGCGATGGCGAACTCAAGGAGCACATCAGCCGCGTCCACGCCGCCAACTACGGTGTTTACGGTGCCCGCAAAG
        TGTGGCTAACCCTGAACCGTGAGGGCATCGAGGTGGCCAGATGCACCGTCGAACGGCTGATGACCAAACTCGGCCTGTCCGGGACCACCCGCGGCAAAGCCC
        GCAGGACCACGATCGCTGATCCGGCCACAGCCCGTCCCGCCGATCTCGTCCAGCGCCGCTTCGGACCACCAGCACCTAACCGGCTGTGGGTAGCAGACCTCA
        CCTATGTGTCGACCTGGGCAGGGTTCGCCTACGTGGCCTTTGTCACCGACGCCTACGCTCGCAGGATCCTGGGCTGGCGGGTCGCTTCCACGATGGCCACCT
        CCATGGTCCTCGACGCGATCGAGCAAGCCATCTGGACCCGCCAACAAGAAGGCGTACTCGACCTGAAAGACGTTATCCACCATACGGATAGGGGATCTCAGT
        ACACATCGATCCGGTTCAGCGAGCGGCTCGCCGAGGCAGGCATCCAACCGTCGGTCGGAGCGGTCGGAAGCTCCTATGACAATGCACTAGCCGAGACGATCA
        ACGGCCTATACAAGACCGAGCTGATCAAACCCGGCAAGCCCTGGCGGTCCATCGAGGATGTCGAGTTGGCCACCGCGCGCTGGGTCGACTGGTTCAACCATC
        GCCGCCTCTACCAGTACTGCGGCGACGTCCCGCCGGTCGAACTCGAGGCTGCCTACTACGCTCAACGCCAGAGACCAGCCGCCGGCTGAGGTCTCAGATCAG
        AGAGTCTCCGGACTCACCGGGGCGGTTCA
        """
        
        self.target = SeqRecord(
            Seq("".join(line.lstrip() for line in target_seq.splitlines())), 
            name='IS6110', id='IS6110', description='')
                             
        self.IRL = SeqRecord(Seq("TGAACCGCCCCGGCATGTCCGGAGACTC"), name='IS6110_IRL', id='IRL', description='')
        self.IRR = SeqRecord(Seq("GAGTCTCCGGACTCACCGGGGCGGTTCA"), name='IS6110_IRR', id='IRR', description='')
        self.IR_length = 28
        self.family = 'IS3'
        self.cluster = 'IS3_168'
        
        self.coordinates = {}
        self.anchors = {}  # insertion IDs as keys, 5' and 3' anchor sequences as values in a tuple
        self.is_seqs = {}  # 
        
        self.IR_blast = []
        self.IR_pairs = []
        self.stats = {
            'properly_paired_IRs': 0,
            'paired_IRs_wrong_internal_seq' : 0,
            'unpaired_IRL': 0,
            'unpaired_IRR': 0
        }
        
    
    def annotate_inverted_repeats(self, path_to_assembly):
        """
        Blast inverted repeats (IRs) of the target IS against the annotated assembly.
        Purpose: identify IRs that are not part of annotated IS and migth produce false positives
        with detettore. Most likely such fragmented IRs are due to poor assembly quality.
        
        """
        # Write IRs to fasta and blast
        with open('query.fa', 'w') as f:
            SeqIO.write(self.IRL, f, 'fasta')
            SeqIO.write(self.IRR, f, 'fasta')
        
        cmd = ['blastn','-query', 'query.fa','-subject', path_to_assembly,
               '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send sstrand pident length']
        
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output = proc.stdout.read()
        
        hits = []
        for line in output.splitlines():
            line = line.decode('utf-8')
            fields = line.split('\t')
            hits.append(fields)

        hits = pandas.DataFrame(hits, columns=['qseqid', 'qstart', 'qend', 'qlen', 'sseqid', 'sstart', 'send', 'sstrand', 'pident', 'length'])
        for var in ['qstart', 'qend', 'qlen', 'sstart', 'send', 'length']:
            hits[var] = hits[var].astype(int)
        hits['pident'] = hits['pident'].astype(float)
        
        # Convert to 0-based index
        hits['sstart'] = hits['sstart'] - 1
        hits['send'] = hits['send'] - 1
            
        self.IR_blast = hits.sort_values(by='sstart', ignore_index=True)
            
        os.remove('query.fa')
        
    
    def evaluate_IS_annotation(self, path_to_assembly, min_id):
        
        assembly = SeqIO.read(path_to_assembly, 'fasta')
        
        i = 0 
        while i < len(self.IR_blast)-1:
    
            row = self.IR_blast.iloc[i]
            row_next = self.IR_blast.iloc[i+1]
                
            ir_start = row['sstart'] if row['sstrand'] == 'plus' else row['send']
            ir_type = row['qseqid']
            
            ir_end_next = row_next['send'] if row_next['sstrand'] == 'plus' else row_next['sstart']
            ir_type_next = row_next['qseqid']
            
            ir_dist = ir_end_next - ir_start + 1
            
            if ir_type != ir_type_next and ((len(self.target)-5) < ir_dist < (len(self.target) + 5)):
                
                strand = '+' if ir_type == 'IRL' else '-'
                seq = assembly[ir_start:ir_end_next+1]
                if strand == '-':
                    seq = seq.reverse_complement()
                dist = sequence_similarity(self.target.seq, seq.seq)      
                
                if dist >= min_id:
                    self.stats['properly_paired_IRs'] += 1
                    self.IR_pairs.append([ir_start, ir_end_next, strand, ir_dist, dist])
                else: 
                    self.stats['paired_IRs_wrong_internal_seq'] += 1    
                i += 2
                
            else:
                self.stats[f'unpaired_{ir_type}'] += 1
                i += 1
                
        
    def get_true_coordinates(self, path_to_isescan_fna):
        """
        Extracts the true coordinates of IS elements from an ISEScan FASTA file.

        This method parses the provided FASTA file, extracts position and strand
        information for each IS element, and calculates sequence similarity with
        the target IS sequence to determine true insertions.
        
        UPDATE: Note overall similarity, but use identity of inverted repeats to
        find "true" insertions. Some of the genomes annotated contain fragmented
        copies that have poor overall similarity, but still produce anchors because
        they have the IRs of the target IS. 
    

        Parameters
        ----------
        path_to_isescan_fna : str
            The path to the FASTA file output by ISEScan, containing IS elements.

        Updates
        -------
        self.coordinates : dict
            Updates the dictionary with sequence record IDs as keys and lists as values.
            Each list contains insertion start position, insertion end position, strand,
            and similarity score rounded to two decimal places.
        self.is_seqs : dict
            Updates the dictionary with sequence record IDs as keys and SeqRecord
            objects as values for the true IS sequences.
        """

        recs = SeqIO.parse(path_to_isescan_fna, 'fasta')
        
        header = ['isbegin', 'isend', 'isstrand', 'target_similarity', 'IRL_similarity', 'IRR_similarity']
        sys.stderr.write('\t'.join(header) + '\n')

        for rec in recs:
            
            ispos = '_'.join(rec.id.split('_')[-3:]).split('_')
            isbegin = int(ispos[0])
            isend = int(ispos[1])
            isstrand = ispos[2] if len(ispos) == 3 else ''
            is_cluster =rec.description.split(' ')[-1]
            
            if is_cluster != self.cluster:
                continue
            
            if isstrand:
                dist = sequence_similarity(self.target.seq, rec.seq)          
                        
            # Infer strand for the elements where no strand was inferred
            else:
                d = sequence_similarity(self.target.seq, rec.seq)
                d_rev = sequence_similarity(self.target.seq, rec.seq.reverse_complement())
                
                if d_rev > d:
                    dist = d_rev
                    isstrand = '-'
                    rec.seq = rec.seq.reverse_complement()
                else:
                    dist = d
                    isstrand = '+'
                
            # Infer similarity of inverted repeats
            IRL_dist = sequence_similarity(self.IRL, rec.seq[:self.IR_length])
            IRR_dist = sequence_similarity(self.IRR, rec.seq[-self.IR_length:])

            self.coordinates[rec.id] = {
                'isbegin': isbegin,
                'isend': isend,
                'isstrand': isstrand,
                'target_similarity': dist,
                'IRL_similarity': IRL_dist,
                'IRR_similarity': IRR_dist
            }
            sys.stderr.write('\t'.join(map(str, list(self.coordinates[rec.id].values()))) + '\n')
 
            self.is_seqs[rec.id] = rec
    
    
    def get_true_anchors(self, assembly_fna, anchor_length, min_id):
        """
        Extracts true anchor sequences from an assembly based on pre-determined IS insertion coordinates.

        Parameters
        ----------
        assembly_fna : str
            The path to the FASTA file containing the assembly sequence.
        anchor_length : int
            The length of the anchor sequences to extract on both sides of the insertion.

        Modifies
        --------
        self.anchors : dict
            Updates the anchors dictionary with record identifiers as keys. Each value is a tuple of two
            SeqRecord objects representing the 5' and 3' anchor sequences, with their orientations and
            descriptions adjusted based on the strand information.
        """

        assembly = SeqIO.read(assembly_fna, 'fasta')
        
        for rec_id in self.coordinates:
            
            if self.coordinates[rec_id]['target_similarity'] < min_id:
                continue
            
            start = self.coordinates[rec_id]['isbegin']
            end = self.coordinates[rec_id]['isend']
            strand = self.coordinates[rec_id]['isstrand']
            
            anchor_front = assembly[(start-anchor_length):(start-1)]  ## gff is 1-based!
            anchor_back = assembly[end:(end+anchor_length)]
            
            # revese complement if necessary to get same anchor orientation as reported by detettore
            if strand == '-':
                anchor_front.seq = anchor_front.seq.reverse_complement()
                anchor_back.seq = anchor_back.seq.reverse_complement()
                
            anchor_front.id = rec_id
            anchor_back.id = rec_id
            anchor_front.description = f'{start-anchor_length}-{start-1}_{strand}'
            anchor_back.description = f'{end}-{end+anchor_length}_{strand}'
                        
            if strand == '+':
                self.anchors[rec_id] = (anchor_front, anchor_back)
            elif strand == '-':
                self.anchors[rec_id] = (anchor_back, anchor_front)
                
                
    def evaluate_detettore_anchors(self, detettore_anchors_tsv, min_id):
    
        """
        Compares the anchor sequences inferred by detettore against the true anchor sequences
        and returns a dictionary with the following keys and values:
        
        copy_number : int
            The number of IS insertions in the assembly.
        true_positives : dict
            A dictionary with '5' and '3' as keys, each containing the number of true positive
            anchor sequences on that side.
        false_negatives : dict
            A dictionary with '5' and '3' as keys, each containing the number of false negative
            anchor sequences on that side.
        false_positives : dict
            A dictionary with '5' and '3' as keys, each containing the number of false positive
            anchor sequences on that side.
        
        Parameters
        ----------
        detettore_anchors_tsv : str
            The path to the TSV file containing the anchor sequences inferred by detettore.
        min_id : float
            The minimum sequence identity required for a true positive match.
        
        Returns
        -------
        dict
            A dictionary with the above keys and values.
        """
        true_anchors = {
            '5' : [self.anchors[x][0] for x in self.anchors],
            '3' : [self.anchors[x][1] for x in self.anchors]
        }
        
        positives = {x : {'5':[],'3':[]} for x in self.anchors}
        negatives = {'5':[],'3':[]}

        with open(detettore_anchors_tsv) as f:
            next(f)
            for line in f:
                fields = line.strip().split('\t')
                side = fields[1]
                rec = SeqRecord(Seq(fields[3]), id=fields[0], name='', description='')
                reclen = len(rec.seq)
                
                # Compare against true anchors! Ignore side information for the moment, just try ever combination to figure out what's going on
                matched = False
                
                for true_anchor in true_anchors[side]:
                    
                    # cut true anchor to same size as inferred
                    if side == '5':
                        true_anchor_seq = true_anchor.seq[-reclen:]
                    else:
                        true_anchor_seq = true_anchor.seq[:reclen]
                                                
                    d = sequence_similarity(rec.seq, true_anchor_seq)
                    
                    # True positive
                    if d > min_id:
                        positives[true_anchor.id][side].append(rec.id)
                        matched = True
                    
                # False positive
                if not matched:
                    negatives[side].append(rec.id)
                    
            outd = {
                'copy_number' : len(positives),
                'true_positives' :  {'5':0,'3':0},
                'false_negatives' : {'5':0,'3':0},
                'false_positives' : {'5':len(negatives['5']),'3':len(negatives['3'])}
            }
            
            for insertion in positives:
                for side in ['5', '3']:
                    if positives[insertion][side]:
                        outd['true_positives'][side] += 1
                    else:
                        outd['false_negatives'][side] += 1
                        
        return outd


#%%
def main():
    
    args = get_args()

    truth = true_insertions()
    truth.get_true_coordinates(args.isescan)
    truth.get_true_anchors(args.assembly, args.anchor_length, args.min_id)
    eval = truth.evaluate_detettore_anchors(args.detettore, args.min_id)
    
    # Check possible problems with IS annotation by blasting inverted repeats
    truth.annotate_inverted_repeats(args.assembly)
    truth.evaluate_IS_annotation(args.assembly, 0.5)
    
    # Write output
    header = [
        'copy_number', 
        'paired_IRs', 
        'paired_IRs_wrong_seq', 
        'unpaired_IRL', 
        'unpaired_IRR',
        'TP_5', 'TP_3', 
        'FN_5', 'FN_3', 
        'FP_5', 'FP_3'
        ]
    
    sys.stdout.write('\t'.join(header)+'\n')
    
    row = [
        eval['copy_number'],
        truth.stats['properly_paired_IRs'],
        truth.stats['paired_IRs_wrong_internal_seq'],
        truth.stats['unpaired_IRL'],
        truth.stats['unpaired_IRR'],
        eval['true_positives']['5'], eval['true_positives']['3'],
        eval['false_negatives']['5'], eval['false_negatives']['3'],
        eval['false_positives']['5'], eval['false_positives']['3']
    ]
    
    sys.stdout.write('\t'.join(map(str,row)) +'\n')
    
if __name__ == '__main__':
    main()