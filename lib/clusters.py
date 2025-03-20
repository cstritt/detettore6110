#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bisect
import os
import pandas
import pysam
import re
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter
from sklearn.linear_model import LinearRegression

from collections import Counter


class AnchorCluster:

    def __init__(self,cluster_nr, side):
        """
        Initialize an AnchorCluster instance with a cluster number and side.

        Parameters
        ----------
        cluster_nr : intipython
        """

        self.cluster_nr = cluster_nr
        self.side = side
        self.cluster_id = f'{side}prime_{cluster_nr}'
        self.reads = []
        self.seqs = []
        self.is_seqs = []  # to store the read parts matching the IS target
        self.ref_coords = ''  # To store pysam read object with .reference_name, .reference_start, .reference_end, .is_reverse, .cigarstring, .mapping_quality
    
    
    def align_anchor_reads(self, temp_dir, cpus):
        
        """
        Align the reads in the AnchorCluster instance with MAFFT.
        
        Parameters
        ----------
        temp_dir : str
            Path to temporary directory.
        args : class
            Input arguments.
        """
        fasta_path = os.path.join(temp_dir, f'{self.cluster_id}.fasta')
        alignment_path = os.path.join(temp_dir, f'{self.cluster_id}.aligned.fasta')
        
        with open(fasta_path, 'w') as fasta_handle:
            SeqIO.write(self.seqs, fasta_handle, 'fasta')
        
        mafft_cmd = [
            'mafft',
            '--thread', str(cpus),
            '--adjustdirection',
            fasta_path
            ]
        
        subprocess.run(mafft_cmd, check=True, 
                    stdout=open(alignment_path, 'w'), 
                    stderr=subprocess.DEVNULL)
        
        os.remove(fasta_path)
        
        
    def summarize_IS_coordinates(self, read_d):
        """ Go through reads and store which positions of the IS are covered
        """
        self.target_cov = {}
        for read_id in self.reads:
            
            read = read_d[read_id]
            for coords in read.target_coordinates:
                start, end = coords
                for i in range(start, end+1):
                    if i not in self.target_cov:
                        self.target_cov[i] = 0
                    self.target_cov[i] += 1
                
        depth = [self.target_cov[i] for i in self.target_cov]
        
        pos = [[i] for i in self.target_cov]  # sort and use index rather than actual position
        
        self.target_lm = lm(depth, pos)
                

    def get_cluster_consensus(self, temp_dir):
        
        """
        Compute the consensus sequence for the cluster based on alignment.
        
        To assess the quality of the cluster, fit an lm(coverage~position). 
        Slope and intercept show if the coverage is increasing or decreasing as 
        expected, or if the cluster is messy and no trend in coverage is there.

        This method reads the alignment from a given file, calculates the
        consensus sequence, and determines the proportion of sites with
        mismatches. It also tracks the sequence depth at the start and end
        of the alignment.

        Parameters
        ----------
        alignment_path : str
            Path to the alignment file in FASTA format.

        Attributes
        ----------
        nr_reads : int
            Number of reads in the alignment.
        aln_len : int
            Length of the alignment.
        consensus : str
            Consensus sequence derived from the alignment.
        depth_start : int
            Depth of sequence coverage at the start of the alignment.
        depth_end : int
            Depth of sequence coverage at the end of the alignment.
        prop_sites_with_mismatches : float
            Proportion of sites with mismatches in the alignment.
        """      
        
        alignment_path = os.path.join(temp_dir, f'{self.cluster_id}.aligned.fasta')

        aln = AlignIO.read(open(alignment_path), "fasta")
        aln_smry = AlignInfo.SummaryInfo(aln)
        
        self.nr_reads = len(aln)
        self.aln_len = aln.get_alignment_length()
        
        consensus = ''
        n_sites_with_mismatches = 0
        
        # For linear model
        depth = []
        position = []
        
        for i in range(self.aln_len):
            col = aln_smry.get_column(i)
            count_missing = col.count('-')
            count_present = self.nr_reads - count_missing
            
            depth.append(count_present)
            position.append([i+1])
            
            # Check on which side the "tail" of the alignment is
            if i == 0:
                self.depth_start = count_present
            if i == (self.aln_len-1):
                self.depth_end = count_present         
            
            # Get consensus base
            if count_present >= 3:
                counter = Counter(col.replace('-', ''))  # don't! Majority can be gap
                base = counter.most_common(1)[0][0]
                
                if len(counter) > 1:
                    n_sites_with_mismatches += 1
            else:
                base = '-'
                
            consensus += base
        
        consensus = consensus.replace('-', '').upper()
        
        self.prop_sites_with_mismatches = round(n_sites_with_mismatches / self.aln_len, 2)

        self.consensus = SeqRecord(
            Seq(consensus),
            id=self.cluster_id,
            name = '',
            description=''
            )
        
        self.anchor_lm = lm(depth, position)


class Clusters:
    def __init__(self, args, temp_dir):
        
        # Input
        self.min_anchor_length = args.min_anchor_len
        self.min_cluster_size = args.min_cluster_size
        self.cpus = args.cpus
        self.temp_dir = temp_dir

        # Output
        self.by_side = {'5': [], '3': []}
        self.cluster_d = {'5': {}, '3': {}}
        
    
    def cluster_anchors(self, reads, seed_len):
        """
        Cluster the anchor sequences in the read dictionary by exact identity.

        Parameters
        ----------
        read_d : dict
            A dictionary of Read objects, keyed by read ID.
        seed_len : int
            The length of the anchor sequence to consider for clustering.

        Attributes
        ----------
        five : dict
            A dictionary of 5' anchor sequences grouped by exact identity.
        three : dict
            A dictionary of 3' anchor sequences grouped by exact identity.
        """
        
        def group_strings_by_value(dictionary):
            """ Group a dictionary by values. Used here to cluster sequences by exact identity
            
            Args:
                dictionary (dict): keys are sequence IDs, values are sequences.

            Returns:
                dict: keys are sequences, values are sequence IDs
            """

            grouped = {}
            for key, values in dictionary.items():
                for value in values:
                    grouped.setdefault(value, []).append(key)
            return grouped
        
        # Clustering, separate for 5' and 3'
        
        identifiers_5 = {}
        identifiers_3 = {}

        for read_id in reads.read_d:
            for rec in reads.read_d[read_id].anchor_seq:
                
                seq = str(rec.seq)
                seq_id = rec.id
                
                if reads.read_d[read_id].side == '5':
                    identifiers_5.setdefault(seq_id, []).append(seq[-seed_len:])
                
                elif reads.read_d[read_id].side == '3':
                    identifiers_3.setdefault(seq_id, []).append(seq[:seed_len])
                    
        self.by_side['5'] = group_strings_by_value(identifiers_5)
        self.by_side['3'] = group_strings_by_value(identifiers_3)
        
        
    def parse_clusters(self, reads):
        """
        Parse clusters of anchor sequences from a dictionary of Read objects: 
            - initiate AnchorCluster
            - align anchor reads and call consensus
            - write consensi to fasta
        
        Addition: also align IS parts and write their consensi to fasta

        Parameters
        ----------
        read_d : dict
            A dictionary of Read objects, keyed by read ID.

        Attributes
        ----------
        parsed : dict
            A dictionary of AnchorCluster objects, keyed by cluster ID and side.
        """
        
        fasta_handle = open(f'{self.temp_dir}/anchor_consensi.fasta', 'w')

        for side in self.by_side:
            
            n = 0
            
            for seq in self.by_side[side]:
                
                # Skip short clusters
                if len(self.by_side[side][seq]) < self.min_cluster_size:
                    continue
                
                ac = AnchorCluster(n, side)
                
                for read in self.by_side[side][seq]: 
                    read_id = read[:-2]
                    read_index = int(read[-1])
                    anchor_rec = reads.read_d[read_id].anchor_seq[read_index]
                    target_rec = reads.read_d[read_id].target_seq[read_index]
                        
                    ac.reads.append(read_id)
                    ac.seqs.append(anchor_rec)
                    ac.is_seqs.append(target_rec)
                                        
                ac.align_anchor_reads(self.temp_dir, self.cpus)
                ac.get_cluster_consensus(self.temp_dir)
                ac.summarize_IS_coordinates(reads.read_d)
                self.cluster_d[side][n] = ac
                    
                if len(ac.consensus) > self.min_anchor_length:
                    SeqIO.write(ac.consensus, fasta_handle, 'fasta')
            
                n += 1
                
        fasta_handle.close()
        
        
    def add_ref_coordinates_to_clusters(self, ref_aligned_anchors):
        """
        Add reference coordinates to each cluster in the cluster dictionary.

        Parameters
        ----------
        ref_aligned_anchors : str
            Path to a BAM file containing aligned anchor sequences.
        """
        pybam = pysam.AlignmentFile(ref_aligned_anchors, "rb")
        
        ref_d = {}

        for read in pybam.fetch():
            ref_d[read.query_name] = read
        pybam.close()
        
        for side in self.cluster_d:
            for cluster_nr in self.cluster_d[side]:
                cluster_id = self.cluster_d[side][cluster_nr].cluster_id
                if cluster_id in ref_d:
                    self.cluster_d[side][cluster_nr].ref_coords = ref_d[cluster_id]

    
    def write_cluster_output(self, args):
    
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
            header += ['target_start', 'target_end', 'prop_sites_with_mismatches']
            
        outhandle.write('\t'.join(header) + '\n')

        for side in self.cluster_d:
    
            for cluster_id in self.cluster_d[side]:
                cl = self.cluster_d[side][cluster_id]
                target_pos = [i for i in cl.target_cov]
                        
                row = [cl.cluster_id, side, cl.nr_reads, str(cl.consensus.seq)]
                
                if args.detailed:
                    if args.reference:
                        try:
                            ref_strand = '-' if cl.ref_coords.is_reverse else '+'
                            row += [
                                cl.ref_coords.reference_name, cl.ref_coords.reference_start, cl.ref_coords.reference_end,
                                ref_strand, cl.ref_coords.cigarstring, cl.ref_coords.mapping_quality
                                ]
                            
                        except AttributeError:
                            row += ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']

                    row += [min(target_pos), max(target_pos),  cl.prop_sites_with_mismatches]

                outhandle.write('\t'.join(map(str, row)) + '\n')
                            
        outhandle.close()


    def write_reference_output(self, overlaps, ref_aligned_reads, args):
    
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
            
        header = ['chromosome', 'position', 'strand', 'TSD', 'support_5', 'support_3','support_ref']
        
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
            header += ['anchor_5', 'anchor_3','mapq_5', 'mapq_3', 'cigar_5', 'cigar_3']
            
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
                
            support_5 = len(self.cluster_d['5'][five_cl_nr].reads)
            support_3 = len(self.cluster_d['3'][three_cl_nr].reads)
            
            # Reference support
            support_ref = get_reference_support(chrom, position, ref_aligned_reads)

            # TSD
            tsd_len = ins[4]
            tsd = self.cluster_d['5'][five_cl_nr].consensus.seq[-tsd_len:]

            # Anchor ID, mapq and cigar
            anchor5_id = ins[0]
            anchor3_id = ins[2]
            mapq_5 = self.cluster_d['5'][five_cl_nr].ref_coords.mapping_quality
            cigar_5 = self.cluster_d['5'][five_cl_nr].ref_coords.cigarstring
            mapq_3 = self.cluster_d['3'][three_cl_nr].ref_coords.mapping_quality
            cigar_3 = self.cluster_d['3'][three_cl_nr].ref_coords.cigarstring

            outline = [chrom, str(position), strand, str(tsd), str(support_5), str(support_3), str(support_ref)]
            
            if args.annot:
                gene, dists_to_gene = gene_overlap(position, annot, chrom_length)
                outline += [gene, dists_to_gene]
                
            if args.detailed:
                outline += [anchor5_id, anchor3_id, str(mapq_5), str(mapq_3), cigar_5, cigar_3]

            outline = map(str, outline)

            outhandle.write('\t'.join(outline) + '\n')
            
        outhandle.close()


def get_reference_support(chromosome, position, bamfile, overlap=20):
    """ Distinguished fixed from non-fixed polymorphisms: find reads that
    overlap the insertion breakpoint
    """
    
    readnr = 0
    
    pybam = pysam.AlignmentFile(bamfile, "rb")
    
    for read in pybam.fetch(chromosome, position, position+1):

        if read.mapq == 0:
            continue

        down = [x for x in range(read.reference_start, read.reference_end) if x < position]
        up = [x for x in range(read.reference_start, read.reference_end) if x > position]

        if len(down) > overlap and len(up) > overlap:
            readnr += 1
            
    pybam.close()
            
    return readnr
        
    



def cd_hit(fasta_path, output_path, min_id = 0.99, ap=False):
    """
    Run cd-hit-est on a fasta file of anchor sequences and return a dictionary
    where keys are cluster numbers and values are lists of read IDs present in
    each cluster.

    NOT USED.
        
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
        '-d', '0',
        '-c', str(min_id),
        '-o', output_path,
        '-sc', '1', 
        '-g', '1'
        ]
        
    if ap:
        cd_hit += ['-ap', '1']  # Doesn't seem to do anything...
        
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

    
def lm(target, features):
    """
    Perform linear regression on depth and positions data.

    This function fits a linear regression model to the given depth and positions
    data, and returns the intercept and slope of the fitted line.

    Parameters
    ----------
    depth : array-like
        The independent variable data (e.g., depth values).
    positions : array-like
        The dependent variable data (e.g., position values).bisect

    Returns
    -------
    intercept : float
        The y-intercept of the regression line.
    slope : float
        The slope of the regression line.
    """

    m = LinearRegression()
    m.fit(features,target)
    intercept = round(float(m.intercept_),2)
    slope = round(float(m.coef_[0]),2)
    return intercept, slope

    

def find_overlaps(ref_aligned_anchors, tsd_len):
    
    """
    Find overlapping reads in a given sorted BAM file.

    Parameters
    ----------
    ref_aligned_anchors : str
        Path to the sorted BAM file containing anchor reads aligned to the reference.
    tsd_len : list
        List of possible target site duplication lengths.

    Returns
    -------
    overlaps : list
        List of tuples containing read IDs, start and end coordinates of overlapping reads,
        and the length of the overlap.
    """
    pybam = pysam.AlignmentFile(ref_aligned_anchors, "rb")
    read_ids = []
    coordinates = []

    for read in pybam.fetch():
        read_ids.append(read.query_name)
        coordinates.append((read.reference_start, read.reference_end))
    pybam.close()

    overlaps = []
    for i in range(len(coordinates)):
        for j in range(i + 1, len(coordinates)):
            start1, end1 = coordinates[i]
            start2, end2 = coordinates[j]
            overlap = end1 - start2
            if overlap in tsd_len:
                if read_ids[i][0] == read_ids[j][0]:  # make sure reads reach into opposite ends of the IS
                    continue
                overlaps.append((read_ids[i], end1, read_ids[j], start2, overlap))
                
    return overlaps
           

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