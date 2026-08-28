#!/usr/bin/env python

""" Create summaries of detettore results:
    - table with all insertions
    - copy number
    - presence-absence matrix
    - frequency spectrum
    - fasta file with all unique anchor sequences
    - strain-specific anchor -> unique anchor map (strain, anchor_id, unique_anchor_id)

To do:
    - apply coverage filter to anchors: remove outliers
    - use cd-hit to identify spurious unique anchors due to indels in the seed
"""
import numpy as np
import os
import pandas
import subprocess
import sys
import tempfile


from Bio import Align
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from scipy.stats import median_abs_deviation

from .readparsing import mapreads
from .clusters import find_overlaps, gene_overlap

#%% Funzioni
def count_mismatches(seq1, seq2, length=20):
    """
    Count the number of mismatched positions between two sequences up to a specified length.
    Used to identify anchors that originate from the same insertion event. 

    Parameters
    ----------
    seq1 : str
        The first sequence to compare.
    seq2 : str
        The second sequence to compare.
    length : int, optional
        The maximum number of positions to compare (default is 20).

    Returns
    -------
    int
        The number of mismatched positions in the compared segment of the sequences.
    """

    mismatches = sum(1 for i in range(min(len(seq1), len(seq2), length)) if seq1[i] != seq2[i])
    return mismatches


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
    mismatches = sum(1 for a, b in zip(best_alignment.target, best_alignment.query) if a != b)
    # Compute the distance as the number of mismatches divided by the length of the alignment
    distance = mismatches / len(best_alignment.target)
    similarity = 1 - distance
    return similarity  


def get_z_score(pandas_column):
    """ 
    Use the z-score to identify coverage outliers.
    
    Conventional thresholds: below or above 2, meaning a value is 
    2 standard deviations below or above the mean.
    
    import numpy as np

    def modified_z_scores(data):
        median = np.median(data)
        mad = median_abs_deviation(data, scale='normal')  # or scale=1.4826 manually
        return 0.6745 * (data - median) / mad

    """
    z_scores = (pandas_column - pandas_column.mean()) / pandas_column.std()
    z_scores = [round(x,2) for x in z_scores]
    
    
    # Modified z-scores for smaller sample sizes
    median = pandas_column.median()
    mad = median_abs_deviation(pandas_column, scale='normal')  # or scale=1.4826 manually
    mod_z_scores = 0.6745 * (pandas_column - median) / mad
    mod_z_scores = [round(x, 2) for x in mod_z_scores]
    
    return z_scores, mod_z_scores


def cd_hit(fasta_path, output_path, min_id = 0.99, ap=False):
    """
    Run cd-hit-est on a fasta file of anchor sequences and return a dictionary
    where keys are cluster numbers and values are lists of read IDs present in
    each cluster.
        
    Parameters
    ----------
    fasta_path : str
        The path to the fasta file to be clustered.
    output_path : str
        The path to the output file (without the .clstr extension).

    Returns
    -------
    clusters : dict
        A dictionary where keys are cluster numbers and values are lists of
        read IDs present in each cluster.
    """
    
    cd_hit = [
        'cd-hit-est',
        '-i', fasta_path,
        '-d', '0',  # length of description, stop at first space
        '-r', '0',  # do only +/+ alignment
        '-c', str(min_id),
        '-o', output_path,
        '-sc', '1',  # sort clusters by size
        '-g', '1'  # cluster into the most similar cluster, rather than the first encountered
        ]
        
    if ap:
        cd_hit += ['-ap', '1']  # Doesn't seem to do anything...
    
    print(' '.join(cd_hit))
    subprocess.run(cd_hit, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Parse output
    clusters = {}
    
    with open(f'{output_path}.clstr') as f:
        for line in f:
        
            if line.startswith('>'):
                cluster_nr = line.strip().split(' ')[-1]
                clusters[cluster_nr] = []
            else:
                read_id = line.strip().split(' ')[1][1:-3]
                clusters[cluster_nr].append(read_id)
    
    return clusters

def mafft(fasta, outfasta):
    """ 
    Align sequences in fasta with mafft. 
    
    Parameters
    ----------
    fasta : str
        Path to fasta file. 
    outfasta :  str
        Path to output file (fasta) for aligned sequences.
    
    Returns
    -------
    
    """
    outhandle = open(outfasta, 'w')
    mafft = ['mafft', '--reorder', '--auto', fasta]
    subprocess.run(mafft, check=True, stdout=outhandle, stderr=subprocess.DEVNULL)
    outhandle.close()
    aln = AlignIO.read(outfasta, 'fasta')
    return aln
    
    
def find_first_last_not_dash(s):
    """
    Given a string, return the first and the last position that are not -.
    Used to identify where the alignment of a nucleotide sequence starts and ends.
    
    Parameters
    ----------
    s : str
        A nucleotide sequence string. 

    Returns
    -------
    first_pos, last_pos : int 
        String index of the first and last non '-' character.
    
    """
    first_pos = -1
    last_pos = -1
    
    # Find first position that is not '-'
    for i in range(len(s)):
        if s[i] != '-':
            first_pos = i
            break
    
    # Find last position that is not '-'
    for i in range(len(s) - 1, -1, -1):
        if s[i] != '-':
            last_pos = i
            break
    
    return first_pos, last_pos
             

def eval_anchor_divergence(reference_insertions, presence_absence_matrix, args):
    """
    For the insertions with identifiable reference position, do we obtain an identical presence-absence vector
    when using anchor sequences as when using reference position and strand?
    
    Output: for each event present in more than one strain, the number of strains,
    the number of strains according to the reference-free approach, the mean number of mismatches
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    
    """
    header = [
        'position', 'strand', 'nr_strains_ref', 'nr_events_5', 'nr_events_3'
        ]
    
    rec_d = {}
    for pos in reference_insertions.unique_insertions:
        if len(reference_insertions.unique_insertions[pos]) > 1:  # Only consider insertions with more than one strain?
            rec_d[pos] = {'5': [], '3': []}
            for strain, anchor5, anchor3 in reference_insertions.unique_insertions[pos]:
                rec_d[pos]['5'].append((strain, anchor5))
                rec_d[pos]['3'].append((strain, anchor3))
    
    ref_anchor_map = {}         
    for pos in rec_d:
        ref_anchor_map[pos] = {}
        for side in ['5', '3']:
            ref_anchor_map[pos][side] = []
            for strain, anchor_id in rec_d[pos][side]:
                anchor = presence_absence_matrix.anchor_d[(strain, anchor_id)]
                ref_anchor_map[pos][side].append(anchor.idx_unique)
        
        nr_strains_ref = len(reference_insertions.unique_insertions[pos])
        nr_events_5 = len(set(ref_anchor_map[pos]['5']))
        nr_events_3 = len(set(ref_anchor_map[pos]['3']))
        
        outline = [pos[0], pos[1], nr_strains_ref, nr_events_5, nr_events_3]
        
        if nr_events_5 == 1:
            nr_strains_5 = len(ref_anchor_map[pos]['5'])
            outline.append(nr_strains_5)
        if nr_events_3 == 1:
            nr_strains_3 = len(ref_anchor_map[pos]['3'])
            outline.append(nr_strains_3)
        print(outline)
    

#%% Classi
    
class ReferenceInsertions:
    def __init__(self):
        
        self.insertions = {}  # strain as key, a list with the detettore output as values
        self.unique_insertions = {}  # with (position, strand) as key, strain and anchor ids as values
        self.anchor_d = {}  # keys: (strain, anchor_id), values: detettore output. Create entries for both the 5' and the 3' anchor!
        
        self.header = [
            'chromosome', 'position', 'strand', 
            'TSD', 
            'support_5', 'support_3',
            'support_ref', 
            'anchor_5', 'anchor_3',
            'mapq_5', 'mapq_3',
            'cigar_5', 'cigar_3'
            ]
        
        self.format = 'default'
        
        self.stats = {
            'total_insertions' : 0,
            'unique_insertion_sites' : set() 
        }
        
    def create_insertion_dict(self, args, samples):
        """
        Walk through the directory given by args.path_detettore_results and open every
        file ending with 'reference_insertions.tsv'. Read the file and append the
        insertions to the dict self.insertions with the strain name as key.
        
        Also fill dictionary with unique insertion sites, with key (position, strand); and the 
        anchor dictionary anchor_d, for later cross-referencing with unique anchors.
        
        
        Parameters
        ----------
        
        args : argparse.Namespace
            Arguments passed through the command line.
        
        Returns
        -------
        
        None
            Fills the following class attributes:
                - 
                - 
        
        """
        
        for root, dirs, files in os.walk(args.path_detettore_results):
            for file in files:
                if file.endswith('reference_insertions.tsv') and not file.startswith('ALL'):
                    
                    strain_pre = file.split('.')
                    strain = '.'.join(strain_pre[:-2])
                    
                    if len(samples) > 0 and strain not in samples:
                        continue
                    
                    self.insertions[strain] = []
                    
                    with open(os.path.join(root, file), 'r') as f:
                        next(f)
                        for line in f:
                            fields = line.strip().split('\t')
                            position = int(fields[1])
                            strand = fields[2]
                            
                            if (position, strand) not in self.unique_insertions:
                                self.unique_insertions[(position, strand)] = []
                            
                            anchor5 = fields[7]
                            anchor3 = fields[8]
                            self.unique_insertions[(position, strand)].append((strain, anchor5, anchor3))
                            
                            self.anchor_d[(strain, anchor5)] = fields
                            self.anchor_d[(strain, anchor3)] = fields
                            self.stats['unique_insertion_sites'].add((position, strand))
                                                            
                            # Check if default or detailed output format
                            if len(fields) != len(self.header):
                                self.format = 'detailed'
                                
                            self.insertions[strain].append(fields)
                            self.stats['total_insertions'] += 1
                                         
        sys.stderr.write(f'{self.stats["total_insertions"]} insertions in {len(self.insertions)} strains.\n')
        sys.stderr.write(f'{len(self.stats["unique_insertion_sites"])} unique insertion sites.')
                            
                            
    def write_ref_insertions(self, args):
        """
        Write the insertions for all strains to a single file.
        
        Parameters
        ----------
        args : class
            Input arguments containing the output file path.
        """
        with open(os.path.join(args.outpath, 'ALL_reference_insertions.tsv'), 'w') as f:
            header = ['strain'] + self.header
            
            if self.format == 'detailed': 
                header += ['gene', 'dist_to_gene']
            
            f.write('\t'.join(header) +'\n')
        
            for strain in self.insertions:
                for insertion in self.insertions[strain]:
                    insertion.insert(0, strain)
                    f.write('\t'.join(insertion) +'\n')
    
    def write_vcf(self, args, samples):
        
        # Create dictionary of unique insertion sites
        
        d_insertions = {}
        for strain in self.insertions:
            for insertion in self.insertions[strain]:
                pos = int(insertion[2])
                
                if pos not in d_insertions:
                    d_insertions[pos] = []
                d_insertions[pos].append(strain)
                
        reference = SeqIO.read(args.reference, 'fasta')
        
        # Write vcf
        with open(os.path.join(args.outpath, 'ALL_reference_insertions.vcf'), 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##reference=../data/MTBC0/MTBC0_v1.1.fasta\n')
            f.write('##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

            header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n'
            f.write(header)

            positions = sorted(d_insertions.keys())
            
            for pos in positions:
                ref_base = reference[pos-1]  
                alt_base = 'TGAACCGCCCCGGCATGTCCGGAGAC'  # doesn't matter, using start of is6110
                genotypes = []
                
                for strain in samples:
                    if strain in d_insertions[pos]:
                        genotypes.append('1')
                    else:
                        genotypes.append('0')
                
                f.write(f'MTBC0\t{pos}\t.\t{ref_base}\t{alt_base}\t.\t.\t.\tGT\t{"\t".join(genotypes)}\n')
                
                
class UniqueAnchor:
    def __init__(self, sequence, id):
        self.id = id
        self.representative = sequence
        self.members = []
        self.zscores = []  # mod_z_scores of the cluster members, accessed in self.write_metadata
        self.cluster = {}  # Does the uanchor cluster with other uanchors? Store cluster id as key and other cluster members as values
        self.context = ''
        

class Anchors:
    """   
    Information to be extracted:
    
        - copy numbers: 
        - presence-absence matrix
        - meta data: 
        - fasta file for each uniqueAnchor
    
    """
    
    def __init__(self):
        
        # filled by self.process_anchors
        self.unique_anchors = {'5': {}, '3': {}}  # Store
        self.strains = []        
        self.anchordf = {}  # store anchor information for each strain (key) as pandas df, with z-scores
        self.anchor_d = {}  # NOT USED? stores the sequence of anchors, with (strain,anchor_id) as key
        self.stats = {
            'nr_unique_sites' : {'5': 0, '3': 0},
            'count_frequencies' : {'5': {}, '3': {}}
        }
        
        # filled by self.check_unique_anchor_redundancy
        self.uanchor_clusters = {}
        self.uanchor_cluster_alignments = {'5': {}, '3': {}}
        
        # filled by self.create_matrix
        self.matrix = ''  # pandas df filled below with create_matrix. 
        
        
    
        self.anchors = {}  # strain as key and a pandas df of the anchors.tsv as value
        
        # filled by self.get_metadata
        self.metadata = ''  # not used?
        self.anchor_metadata = {}  # reference position and genomic context, if there is, retrieved from ReferenceInsertions
 
                    
                    
    def process_anchors(self, seed_len, max_mismatches, args, samples, min_num_reads=5, min_z_score=-2, subsample=False): 
        """ 
        Loop through anchor.tsv files, calculate z-scores for the read support, combine all anchors into single file,
        find unique anchors by matching the first seed_len bases of the anchors.
        
        Keep the longest anchor sequence as the reference for a UniqueAnchor.
        
        Add filter flag: min_num_reads, min_z_score.
        
        Ignore the reference information, which is also present in the ReferenceInsertions.        
        
        Parameters
        ----------
        seed_len : int
            If anchors match over the first seed_len bases, they reflect the same unique insertion event.
            
        max_mismatches : int
            Number of mismatches allowed for a match.
            
        args : 
        
        samples : list
            Only consider samples listed here, rather than all present in the results directory. 
        
        outfile : str
            Path to file into which all anchors are combined.

        subsample : int
            Only consider n anchor files, to reduce waiting time for development.
            
        Returns
        -------
        
        Creates outfile with all anchors and z-score values.
        Fills self.unique_anchors with members of UniqueAnchor, separately for 5' and 3' sites
        
        """
        outhandle = open(os.path.join(args.outpath, 'ALL_anchors.tsv'), 'w')
        cnt = 0  # for subsampling a large data set
        
        drop_cols = ['ref', 'ref_start', 'ref_end', 'ref_strand', 'ref_cigar', 'ref_mapq']
        
        for root, dirs, files in os.walk(args.path_detettore_results):
            for file in files:
                
                if file.endswith('.anchors.tsv'):
                    
                    strain_pre = file.split('.')
                    strain = '.'.join(strain_pre[:-2])
                    
                    if len(samples) > 0 and strain not in samples:
                        continue
                    
                    self.anchor_d[strain] = {}
                    
                    if subsample and cnt == subsample:
                        return
                        
                    sys.stdout.write(f'Processing {strain} - {cnt}\n')
                    try:
                        anchor_df = pandas.read_csv(os.path.join(root, file), sep='\t')
                    except:
                        continue
                    
                    #nrow = anchor_df.shape[0]
                    #strain_col = pandas.Series(nrow*[strain], name='strain')
                    #anchor_df = pandas.concat([strain_col.to_frame().T, anchor_df], axis=1)
                    anchor_df['strain'] = strain
                    self.strains.append(strain)
                    
                    for column in drop_cols:
                        if column in anchor_df.columns:
                            anchor_df = anchor_df.drop(column, axis=1)
                    
                    # Add z-scores
                    z_scores, mod_z_scores = get_z_score(anchor_df['num_reads'])
                    anchor_df['z_score'] = z_scores
                    anchor_df['mod_z_score'] = mod_z_scores
                    
                    # Add filter flag: PASS and FAILED
                    anchor_df['filter'] = 'PASS'
                    filter_failed = (anchor_df['num_reads'] < min_num_reads) | (anchor_df['mod_z_score'] < min_z_score)
                    #anchor_df['filter'][filter_failed] = 'FAILED'
                    anchor_df.loc[filter_failed, 'filter'] = 'FAILED'
                    
                    # Append all anchors to file
                    is_first = True if cnt == 0 else False
                    anchor_df.to_csv(outhandle, mode='a', sep='\t', header= is_first, index=False)
                    
                    # Remove filtered. They will be in the ALL_anchors file, but not in the matrix
                    anchor_df = anchor_df[anchor_df['filter'] != 'FAILED']

                    # Get information for copy number estimation
                    self.anchordf[strain] = anchor_df
                
                    # Evaluate anchor uniqueness    
                    self.check_uniqueness(strain, anchor_df, seed_len, max_mismatches)
                    
                    cnt += 1
                            
        outhandle.close()
                                  
                                  
    def check_uniqueness(self, strain, anchor_df, seed_len, max_mismatches):
        
        for i, anchor in anchor_df.iterrows():
                        
            anchor_side = str(anchor['side'])
            anchor_seq = anchor['consensus']
            mod_zscore = anchor['mod_z_score']
                                    
            # Compare to existing unique anchors
            best_match = ('', float('inf'))
        
            for ua_id, uanchor in self.unique_anchors[anchor_side].items():
            
                if anchor_side == '5': 
                    mismatches = count_mismatches(anchor_seq[-seed_len:], uanchor.representative[-seed_len:])
                else:
                    mismatches = count_mismatches(anchor_seq[:seed_len], uanchor.representative[:seed_len])
                
                if mismatches < best_match[1]:
                    best_match = (ua_id, mismatches)
                    
            best_id, best_n_mismatches = best_match
                
            # New unique anchor
            if best_n_mismatches > max_mismatches:
                
                uanchor_id = f'{anchor_side}.ua{len(self.unique_anchors[anchor_side])}'
                new_unique = UniqueAnchor(anchor_seq, uanchor_id)
                new_unique.members.append((strain, anchor['anchor_id']))
                new_unique.zscores.append(mod_zscore)
                self.unique_anchors[anchor_side][uanchor_id] = new_unique
                #anchor['idx_unique'] = len(self.unique_anchors[anchor_side]) - 1  #?
                
            # Add to existing 
            else:
                self.unique_anchors[anchor_side][best_id].members.append((strain, anchor['anchor_id']))
                self.unique_anchors[anchor_side][best_id].zscores.append(mod_zscore)
                if len(anchor_seq) > len(self.unique_anchors[anchor_side][best_id].representative):
                    self.unique_anchors[anchor_side][best_id].representative = anchor_seq
                #anchor['idx_unique'] = best_id
                
            self.anchor_d[strain][anchor['anchor_id']] = anchor
    
    
    def cluster_unique_anchors(self, min_id=0.9):
        """ Multiple unique anchors could be created if there are indels within the seed considered for matching anchors.         
        To fix this, use cd-hit to cluster unique anchors and find anchors with the same starting positions.       
        
        Add cluster_id and cluster members to each unique anchor
        
        Parameters
        ----------
        min_id : int
        
        
        Returns
        -------
        None.
            Fills self.uanchor_cluster_alignments
        
        """
        
        tmpdir = tempfile.mkdtemp()
        
        seqd = {}
        alnd = {}
        
        for side in self.unique_anchors:
            
            # Write seqs to fasta
            outfasta = os.path.join(tmpdir, f'unique_anchors.{side}prime.fasta')
            outhandle = open(outfasta, 'w') 
            
            for ua_id, uanchor in self.unique_anchors[side].items():
                
                rec = SeqRecord(Seq(uanchor.representative), name=ua_id, id=ua_id, description='')
                SeqIO.write(rec, outhandle, 'fasta')
                seqd[ua_id] = rec
                
            outhandle.close()
            
            # Cluster them with cd-hit
            outcdhit = os.path.join(tmpdir, f'unique_anchors.{side}prime.cdhit')
            self.uanchor_clusters[side] = cd_hit(outfasta, outcdhit, min_id = min_id)
            
            # Add cluster information to unique_anchors
            for cluster_nr in self.uanchor_clusters[side]:
                members = self.uanchor_clusters[side][cluster_nr]
                for ua_id in members:
                    self.unique_anchors[side][ua_id].cluster = {cluster_nr : [x for x in members if x != ua_id]}
            
            # Align cluster sequences
            for cluster_nr in self.uanchor_clusters[side]:
                if len(self.uanchor_clusters[side][cluster_nr]) < 2:
                    continue
                
                # Write cluster fasta and align
                cluster_recs = [seqd[x] for x in self.uanchor_clusters[side][cluster_nr]]
                cluster_fasta = os.path.join(tmpdir, f'{side}.cluster_{cluster_nr}.fasta')
                cluster_fasta_outhandle = open(cluster_fasta, 'w')
                SeqIO.write(cluster_recs, cluster_fasta_outhandle, 'fasta')
                cluster_fasta_outhandle.close()
                
                # Keep fastas for visual inspection
                #aln = mafft(cluster_fasta, f'/scicore/home/gagneux/stritt0001/TB/random_stuff/anchors/{side}.cluster_{cluster_nr}.aligned.fasta')
                aln = mafft(cluster_fasta, os.path.join(tmpdir, f'{side}.cluster_{cluster_nr}.aligned.fasta'))
                self.uanchor_cluster_alignments[side][cluster_nr] =  aln

    
    def create_matrix(self, side, args, to_file=False):
        """ Create matrix with unique insertion sites as rows and samples as columns.
        Also create metadata file, where insertion IDs are linked to the representatitve anchor sequence,
        the number of strains with presence, and the best blast hit.
        
        Parameters
        ----------
        side : str
        args : argparse.Namespace
        to_file : bool
        
        Returns
        -------
        None
            Fills self.matrix with a pandas DataFrame
        
        """
        self.matrix = pandas.DataFrame(index=self.strains, columns=[x for x in self.unique_anchors[side]])

        # Add samples as columns
        for ua_id, ua in self.unique_anchors[side].items():
            strains_present = [x[0] for x in ua.members]
            self.matrix[ua_id] = [1 if strain in strains_present else 0 for strain in self.strains]

        # Save matrix to TSV file
        if to_file:
            self.matrix.to_csv(os.path.join(args.outpath, f'ALL_presence-absence.{side}prime.tsv'), sep='\t')

        
    def get_frequency_spectrum(self, args, to_file=False):
        col_sums = self.matrix.sum(axis=0)
        self.frequency_spectrum = col_sums.value_counts().sort_index()
        if to_file:
            self.frequency_spectrum.to_csv(os.path.join(args.outpath, 'ALL_frequencies.tsv'), sep='\t', header=None)
   
   
    def characterize_anchors(self, args, reference_insertions):
        """ Blast representative anchor sequences against cds and intergenic regions
        """
        
        anchor_representatives = [ua.representative for ua in self.unique_anchors]
        pass
    
    def estimate_copy_numbers(self, outpath):
        """
        Estimate IS copy numbers. 
        CN is the smaller value of the number of 5' vs 3' anchors.
        For CN_zfilt observations with a modified z-score smaller than -2 are ignored.
        
        anchor_df[['side', 'num_reads', 'z_score', 'mod_z_score']]
        
        Parameters
        ----------
        outpath : str
            Path to output directory.
        min_n_reads : int
            Minimum number of supporting reads to be considered for copy number estimation.
        min_z_score : 
            Minimum adjusted z-score 
            
        Returns
        -------
        
        """
        
        header = ['strain', 'CN', 'nr_5_anchors', 'nr_3_anchors']

        with open(os.path.join(outpath, 'ALL_copy_numbers.tsv'), 'w') as f:
    
            f.write('\t'.join(header) +'\n')
            for strain in self.anchordf: 
                
                anchor_df = self.anchordf[strain]
                anchor_df_5 = anchor_df[anchor_df['side'] == 5]
                anchor_df_3 = anchor_df[anchor_df['side'] == 3]
                
                n5 = anchor_df_5.shape[0]
                n3 = anchor_df_3.shape[0]
                cn = min(n5, n3)
                
                f.write(f'{strain}\t{cn}\t{n5}\t{n3}\n')
                
                
    def write_metadata(self, args, reference_insertions):
        """ Write meta data for the unique anchors. Also write fasta file for each unique anchor, 
        containing the sequences of all its members.      
        
        Also match 5' and 3' anchors if they map to a reference.
        
        + Write a separate file mapping anchor IDs to reference insertions
        
        def rms_zscore_quality(z_scores):
            Calculate RMS of z-scores as an overall quality metric.
            return np.sqrt(np.mean(np.square(z_scores)))
        
        Parameters
        ----------
        
        Returns
        -------
        
        """
        header = [
            'site_id', 
            'side', 
            'ref_position', 
            'ref_strand',
            'context', 
            'ref_position_alt',
            'nr_strains', 
            'cluster',
            'rms_zscores',
            'mean_zscore'
            #'members', 'zscores'
            ]

        with open(os.path.join(args.outpath, 'ALL_presence-absence.metadata.tsv'), 'w') as f:
            
            f.write('\t'.join(header) +'\n')
            
            recs = []

            for side in ['5', '3']:
                
                for site_id in self.unique_anchors[side]:
                    
                    ua = self.unique_anchors[side][site_id]
                    repres = ua.representative
                    nr_strains = len(ua.members)
                    #members = ','.join([f'{x[0]}_{x[1]}' for x in ua.members])
                    
                    #zscores = ','.join(map(str,ua.zscores))
                    rms_zscores = np.sqrt(np.mean(np.square(ua.zscores)))
                    mean_zscore = sum(ua.zscores) / nr_strains
                    
                    # Clustering with other unique anchors?
                    clusters = [cluster_nr for cluster_nr in ua.cluster]
                    clusters_string = ','.join(map(str,clusters))

                    # Get majority reference position, if there is one
                    positions = {}

                    for anchor in ua.members:
                        if anchor in reference_insertions.anchor_d:
                            fields = reference_insertions.anchor_d[anchor]
                            
                            position = fields[2]
                            strand = fields[3]
                            context = fields[-2]
                            
                            positions[(position, strand, context)] = positions.get((position, strand, context), 0) + 1
                            
                    # Find entry with highest count
                    if positions:
                        position, strand, context = max(positions, key=positions.get)
                        positions_alt = ','.join([p for p,s,c in positions.keys() if p != position])
         
                    else:
                        position = 'NA'
                        context = 'NA'
                        strand = 'NA'
                        positions_alt = 'NA'
                    
                    # Write output
                    outline = [
                        site_id, 
                        side, 
                        position, 
                        strand,
                        context,
                        positions_alt, 
                        nr_strains,
                        clusters_string,
                        round(rms_zscores,2),
                        round(mean_zscore,2)
                        ]
                    
                    f.write('\t'.join(map(str, outline)) + '\n')
                    
                    # Fasta entry
                    seqname = f'{side}_{site_id}'
                    rec = SeqRecord(Seq(repres), id = seqname, name = seqname, description = f'nstrains={nr_strains};position={position};context={context}')
                    recs.append(rec)
    
            SeqIO.write(recs, os.path.join(args.outpath, f'ALL_unique_anchors.fasta'), 'fasta')
                
                
    def write_ua_sequences(self, args):
        """ Write a fasta for each unique anchor, containing all member anchors. 
        Useful to inspect how complex, repetitive regions might affect the
        identification of unique anchors and the inference of reference positions. 
        
        Also create combined fasta with all sequences. Used for mapping against 
        reference (map_anchors_vs_reference function below).
        
        Parameters
        ----------
        
        Returns
        -------
        
        """
        try:
            os.mkdir(os.path.join(args.outpath, 'unique_anchors'))
        except FileExistsError:
            pass
                                  
        for side in ['5', '3']:
            
            for site_id in self.unique_anchors[side]:
                
                recs = []                   
                ua = self.unique_anchors[side][site_id] 
                for member in ua.members:
                    strain, anchor_id = member
                    seq = self.anchor_d[strain][anchor_id]['consensus']
                    rec = SeqRecord(
                        Seq(seq), 
                        id = f'{strain}_{anchor_id}', 
                        name = f'{strain}_{anchor_id}', 
                        description = '')
                    recs.append(rec)
                    
                SeqIO.write(recs, os.path.join(args.outpath, 'unique_anchors', f'{site_id}.fasta'), 'fasta')
        
                    
        
    def write_anchor_map(self, args):
        """ Write a file mapping strain-specific anchors to unique anchor IDs.
        """
        with open(os.path.join(args.outpath, 'anchor_map.tsv'), 'w') as f:
            for side in ['5', '3']:
                for site_id in self.unique_anchors[side]:
                    ua = self.unique_anchors[side][site_id]
                    for member in ua.members:
                        strain, anchor_id = member
                        f.write(f'{strain}\t{anchor_id}\t{site_id}\t{side}\n')
                        
