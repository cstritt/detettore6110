#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Classes and functions for detettoreTB

"""

import subprocess
import pysam
import re
import numpy
import os
import sys
import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter



#%% classi


class mise_en_place:
    """ 
    Object for easy access to program settings and input files

    """

    def __init__(self, args):

        # Files
        self.FASTQ = [os.path.abspath(x) for x in args.FASTQ]
        self.REF = os.path.abspath(args.REF)
        self.TARGETS = os.path.abspath(args.TARGETS)
        self.ANNOT = os.path.abspath(args.ANNOT)
        self.OUTPREF = args.OUTPREF
        self.CPUS = args.CPUS
        
        # Temporary folder
        self.tmp = self.OUTPREF + '_tmp'

        # Thresholds
        self.MAPQ = args.MAPQ
        self.MIN_LEN = args.MIN_LEN
        self.MAX_LR_DIST = args.MAX_LR_DIST
        
        # Read info
        self.readlength = int()

        # IS target library
        self.telib = {
            seq_record.id:
                seq_record.seq for seq_record in SeqIO.parse(self.TARGETS, "fasta")}
            
        # Reference chromosomes
        self.chromosomes = [seq_record.id for seq_record in SeqIO.parse(self.REF, 'fasta')]

        # Load reference TE annotation (ISEScan gff format assumed)
        self.annotation = []

        # Check file format
        formt = self.ANNOT.split('.')[-1]
        if not formt.startswith('gff'):
            sys.exit('Check annotation format.')

        with open(self.ANNOT) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                
                # Skip TIRs present in isescan gff files                    
                if fields[2] == 'terminal_inverted_repeat':
                    continue
                
                entry = {
                    # GFF uses 1-based indexing!
                    'chromosome' : fields[0],
                    'start' : int(fields[3]),
                    'end' : int(fields[4]),
                    'strand' : fields[6],
                    'id' : re.search(r'ID=(.*?)[,;]', fields[-1]).group(1),
                    'family' : re.search(r'family=(.*?)[,;]', fields[-1]).group(1),
                    'cluster' : re.search(r'cluster=(.*?)$', fields[-1]).group(1)
                    }

                self.annotation.append(entry)
                
    def mask_ref(self, tmp):
        """ 
        Create reference in which IS annotated in the reference are masked.
        Uses bedtools maskfasta. 
        
        """
        
        with open(tmp + '/refIS.bed', 'w') as f:
            
            for IS in self.annotation:
                
                if IS['cluster'] == 'IS3_168':
                
                    bedline = '\t'.join([IS['chromosome'], str(IS['start']), str(IS['end'])])
                    f.write(bedline + '\n')


        subprocess.run(
            ('bedtools', 'maskfasta', '-fi', self.REF, '-bed', tmp + '/refIS.bed', '-fo', tmp + '/ref_masked.fasta')
            )

            


class read_cluster:
    """
    Container for read clusters supporting non-reference TE insertion. Here
    mapping results are combined with the anchor reads.
    
    """

    def __init__(self, reads):

        self.chromosome = str()
        self.position = int()
        self.IQR = int()

        # Set with names of clustering reads
        self.reads = reads
        self.te_hits = {}
        self.breakpoint = [ [] , [] ]
        self.region_start = float('inf')
        self.region_end = -float('inf')

        # Mapping quality of reads bridging the supposed TE insertion site,
        # thus providing evidence against an insertion. Used to calculate genotype quality
        self.ref_support = {}



class TIP:

    def __init__(self):
    
            self.chromosome = str()
            self.position = int()
    
            # Set with names of clustering reads
            self.reads = []
            self.te_hits = {}
            self.breakpoint = [ [] , [] ]
            self.region_start = float('inf')
            self.region_end = -float('inf')
    
            # Mapping quality of reads bridging the supposed TE insertion site,
            # thus providing evidence against an insertion. Used to calculate genotype quality
            #self.ref_support = {}
            
    def evaluate_reads(self, reads, hits):
        
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
                    'aligned_positions' : set(),
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
            self.te_hits[target]['aligned_positions'].update(
                range(hit['target_start'], hit['target_end'])
                )


#%% funzioni

def mapreads(FASTQ, REF, OUTPREF, CPUS, OUTFMT, tmp, k=15, m=40):   
    
    """ Map Illumina reads against a reference. 
    Minimap2 is used with -Y (use softclipping for supplementary alignments)
    
    Parameters:
        OUTFMT: output format, either sam or paf

    
    """
    
    CPUS = str(CPUS)    
    
    
    if OUTFMT == 'sam':
        
        minimap_cmd =  [
            'minimap2', '-a', '-x', 'sr', '-Y','-t', CPUS, '-k', str(k), '-m', str(m), REF
        ]
    
    elif OUTFMT == 'paf':
        
        minimap_cmd =  [
            'minimap2', '-x', 'sr', '-c', '-o', tmp + '/' + OUTPREF + '.paf', '-Y','-t', CPUS, '-k', str(k), '-m', str(m), REF
        ]

    else:
        return 'Error: specify output format'
    
    minimap_cmd += FASTQ
    #print(' '.join(minimap_cmd))
    
    
    if OUTFMT == 'sam': # Get sorted and indexed bam file
        
        minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        
        subprocess.check_output(
            ('samtools', 'view', '-@', CPUS, '-hu', '-o', OUTPREF + '.raw.bam'),
            stdin=minimap.stdout, stderr=subprocess.DEVNULL)
        
        minimap.wait()
        
        subprocess.check_call(
            ('samtools', 'sort', '-@', CPUS, OUTPREF + '.raw.bam', '-o', OUTPREF + '.bam'), stderr=subprocess.DEVNULL
            )
    
        subprocess.check_call(
            ("samtools", "index", OUTPREF + '.bam')
            )
        
        subprocess.check_output(
            ('rm', OUTPREF + '.raw.bam')
            )
        
    elif OUTFMT == 'paf':
        
        subprocess.run(minimap_cmd, stderr=subprocess.DEVNULL)
    
    
def getsplitreads(BAM, tmp, MIN_LEN=10, MAPQ=1):
    """ Extract splitreads and write split parts to fasta file
    
    Issue: many reads encountered twice...
    -> use read.tags, SA
    

    SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ Other canonical alignments in a chimeric alignment, 
    formatted as a semicolon-delimited list. Each element in the list represents a part of the chimeric 
    alignment. Conventionally, at a supplementary line, the first element points to the primary line. Strand is
either ‘+’ or ‘-’, indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based
coordinate.

    NM: number of mismatches
    

    TO ADD:
        Also write the sequences that map to the genome to a fasta. 
        These can then be clustered with cd-hit to obtain a more accurate estinante
        if the IS6110 copy number. 




    """

    pybam = pysam.AlignmentFile(BAM, "rb")

    splitreads = {}
    
    splitseqs = []

    for read in pybam.fetch():
        
        cigar = read.cigartuples
        
        if len(cigar) > 1:
    
            first, last = cigar[0], cigar[-1]
            
            # what if both ends are clipped? like this, first should be called
            # more often
            if first[0] in (4,5):
                cliplen = first[1]

            elif last[0] in(4,5):
                cliplen = last[1]
                
            else:
                continue
            
            if cliplen < MIN_LEN:
                continue
            
            if read.mapping_quality < MAPQ:
                 continue
             
            # if read.is_supplementary:
            #     continue

            # clseq = read.query_sequence
            # if clseq.count('N')/float(len(clseq)) > 0.1:
            #     continue

            # Write clipped to fasta
            name = read.query_name
            if read.is_read1:
                name = name +'/1'
            elif read.is_read2:
                name = name + '/2'
            
            seqrec = SeqRecord(
                Seq(clip_seq(read.query_sequence, cigar)), 
                id = name, 
                name="", 
                description="")
            
            if name not in splitreads:
                splitreads[name] = []
             
            splitreads[name].append(read)
            
            splitseqs.append(seqrec)
        
    pybam.close()
    
    SeqIO.write(splitseqs, tmp + '/softclipped.fasta', 'fasta')

    return splitreads



def ISmap(queries, targets, outpref, min_aln_len, k, w):

    """ Map discordant reads and split reads against TE consensus library. Mapq is the crucial
    output here, as it is later used to calculate genotype likelihoods.
    
    k: Minimizer k-mer length
    w: Minimizer window size [2/3 of k-mer length]. A minimizer is the smallest k-mer in a window of w consecutive k-mers.
    
    If TEs contain internal repeats, e.g LTR retrotransposons, this can result in multimappers and
    a mapq of 0. To correct for this, mapq is recalculated, setting s2 (f2 in Li 2018) to 0.

    mapq = 40 * (1-s2/s1) * min(1,m/10) * log(s1)

    Li 2018: "Minimap2: Pairwise alignment for nucleotide sequences"
    """

    cmd = [
        'minimap2', 
        '-xsr', 
        '--secondary=yes', 
        '-k', str(k), 
        '-w', str(w),
        targets, 
        queries
        ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output_raw = proc.stdout.read().splitlines()

    outfile = open('%s.k%i.paf' % (outpref, k), 'w')

    minimapd = {}

    for line in output_raw:

        line = line.decode('utf-8')
        outfile.write(line + '\n')
          
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

        # Secondary alignments: if they align to the same element, recalibrate mapping quality.
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


def clip_seq(seq, cigar):
    """
    Return clipped part of sequence, as indicated in cigar. Assumes that there
    is only one clipped part at the beginning or end of the read. Cigar input
    parameter is a pysam cigartuples object.
    """
    if cigar[0][0] == 4:
        end = cigar[0][1]
        clipseq = seq[:end]
    elif cigar[-1][0] == 4:
        strt = sum([x[1] for x in cigar[:-1]])
        end = strt + cigar[-1][1]
        clipseq = seq[strt:end]
    return clipseq


def remove_outliers(lista):
    # outliers defined as in R's boxplots
    q_1 = numpy.percentile(lista, 25)
    q_3 = numpy.percentile(lista, 75)

    boxlength = q_3 - q_1

    lower = max(min(lista), q_1 - 1.5*boxlength)
    upper = min(max(lista), q_3 + 1.5*boxlength)

    lista_filt = [x for x in lista if x >= lower and x <= upper]
    return lista_filt


def consensus_from_bam(region, bamfile, filters):

    """ Create a pileup file for a region and extract consensus sequence
    Quality filtering:
            https://gist.github.com/amblina/5ed9a61ce74ad90668f4e29d62e7eb79
    """

    pybam = pysam.AlignmentFile(bamfile, "rb")

    chrmsm, strt, end = region[0], region[1]-1, region[2]
    baseq, mapq = filters

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

            elif read.mapping_quality < mapq:
                crap_reads.add(read.query_name)
                continue

            elif read.query_qualities[pileupread.query_position] < baseq:
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


def TE_hit_summary(reads, hits):
    """ For a set of clustering reads, create a summary of their
    mapping against the TE library. 
    
    PARAMETERS:
        reads: a list of reads
        hits: dictionary output of strumenti.ISmap
    
    RETURNS:
        A summary dictionary and the name of the IS with the highest score
        (cummulative mapping quality)
    
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
    
    PARAMETERS:
        splitreads: dictionary {readname:[reads], ...}, as output by 
            strumenti.getsplitread()
        hits: dicitonary {readname: {te hit dictionary}}, as output by 
            strumenti.ISmap
    
    
    RETURNS:
        dictionary {(position, side):[reads]}
    
    
    """
    
    # Cluster split reads
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


    
    
    
def merge_clusters(split_positions, hits, TSD_MAX = 10, MAX_DIST = 10000):
    """ Go through clusters of split reads and merge them when they point to 
    the same insertion event.
    
    In case of precise insertions, L and R clusters will overlap. 
    
    Imprecise insertion sites arise when SVs have messed up the insertion signature.
    Often this results in L and R clusters that are separated, e.g.:
    
    MTBC0	1993886	-		L	16	0	16	1258	1353
    MTBC0	1994193	-		R	0	8	8	1	81
    
    
    Parameters
    ----------
    split_positions : dict
    
    hits : dict
    
    TSD_MAX : int
    MAX_DIST : int
    
        
    Returns
    -------
    
    
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
        dist = abs(next_pos[0] - pos[0])
        
        
        # Precise insertion with TSD
        if pos[1] == "right" and next_pos[1] == "left" and dist <= TSD_MAX:
            reads = split_positions[pos] + split_positions[next_pos]
            cat = "LR"
            skip = True
            
        
        # Imprecise: L and R clusters with large or no overlap
        # Assume same insertion event if same strand and opposite ends of the IS
        # are covered
        
        elif (pos[1] == "right" and next_pos[1] == "left" and dist > TSD_MAX) or \
            (pos[1] == "left" and next_pos[1] == "right" and dist < MAX_DIST):
                
            
            if is_same_event(
                    (split_positions[pos], split_positions[next_pos]),
                    hits
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
        
        cl = TIP()
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
        DESCRIPTION.
    hits : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

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
    


def get_partially_mapping(PAF, FASTQ, args, min_anchor_len=20):
    """ 
    Extract reads partially mapping against IS, given minimap2 paf output.    

    Parameters
    ----------
    PAF : str
        Path to PAF output of minimap, containing reads that map to IS.
    FASTQ : str
        Fastq file with reads that were mapped.
        
    min_anchor_len : int
        Minimum length of the anchor/clipped read part
        

    Returns
    -------
    Dictionary with read IDs as keys. Values: the side of the IS to which the
    read maps, the start and the end coordinates of the anchor/clipped read part

    """
    
    partial = {}

    with open(PAF) as f:
        
        for line in f:
            
            fields = line.strip().split('\t')
            
            query_len = int(fields[1])
            aln_len = int(fields[10])
            
            nonmap_len = query_len - aln_len
            
            if nonmap_len > min_anchor_len:
                
                readname = fields[0]
                
                # To which side of the IS does the read map?
                target_start = int(fields[7])
                target_end = int(fields[8])
                target_len = int(fields[6])
                
                query_start = int(fields[2])
                query_end = int(fields[3])
                
                if target_start == 0:
                    side = '5'

                elif target_end == target_len:
                    side = '3'
                    
                else:
                    continue
                
                if query_start == 0:
                    anchor_start = query_end
                    anchor_end = query_len
                
                elif query_end == query_len:
                    anchor_start = 0
                    anchor_end = query_start       
                    
                else:
                    continue
                
                partial[readname] = (side, anchor_start, anchor_end)
                        
    return partial
                    

def subset_fastq(partially_mapping, FASTQ, params):
    """
    Extract the partially mapping reads from the original fastq

    Parameters
    ----------
    partially_mapping : dict
        Output of the get_partially_mapping function.
    FASTQ : list
        Original FASTQ file(s).
    args : class
        Input arguments.

    Returns
    -------
    One (SE) or two (PE) fastq files containing the IS-mapping reads,
    two fasta files containing the anchors of the 5' and the 3' sides.

    """
    
    fastq_out = []

    fasta_out = {
        '5' : [],
        '3' : []
        }
    
    anchor_d = {}
    
    for f in FASTQ: 
        
        with gzip.open(f, "rt") as fastq_handle:
        
            for read in SeqIO.parse(fastq_handle, "fastq"):
                
                if read.id in partially_mapping:
                    
                    fastq_out.append(read)
                    
                    side = partially_mapping[read.id][0]
                    
                    anchor_start = partially_mapping[read.id][1]
                    anchor_end = partially_mapping[read.id][2]
                    anchor_seq = read.seq[anchor_start:anchor_end]
                    
                    fasta_rec = SeqRecord(
                        anchor_seq,
                        id=read.id,
                        name = read.id,
                        description=side + '_' + str(anchor_start) + '-' + str(anchor_end)
                        )

                    fasta_out[side].append(fasta_rec)    
                    anchor_d[read.id] = anchor_seq
             
    with gzip.open(params.tmp + '/partially_mapping.fastq.gz', 'wt') as fastq_handle:
        SeqIO.write(fastq_out, fastq_handle, 'fastq')
    
    for side in fasta_out:
        with open(params.tmp + '/anchors.' + side + '.fasta', 'w') as fasta_handle:
            SeqIO.write(fasta_out[side], fasta_handle, 'fasta')
            
    return anchor_d
        

def remove_known_insertions(clusters, params, imprec = 5):

    """ Develop this part by using RW-TB008 as a reference, so each hit should be 
    present.

    Conserved IS only appear in the output when mapping against masked assembly;
    when mapping agains unmasked, conserved insertions don't produce split reads.


    """    

    remove = []
    
    for i,c in enumerate(clusters):
        
        pos = c[0]
        side = c[1]
        
        for element in params.annotation:
            
            if element['cluster'] != 'IS3_168':
                continue
            
            if side == "L":
                
                if element['start'] - imprec <= pos <= element['start'] + imprec:                
                    print('next to annotated:', pos, side)
                    remove.append(i)
                    
            elif side == "R":
                
                if element['end'] - imprec <= pos <= element['end'] + imprec:
                    print('next to annotated:', pos, side)
                    remove.append(i)
    
    return remove
                

def write_output(params, clusters, tmp, to_file=False, to_list=False):
    """ Writes to stdout by default. 
    
    Similar to ISMapper output.
    
    Add: 
        - mapQ_left, mapQ_right
        - dist L and R support
        
        
    Default filters??:
        - only insertions with support from both sides
        - maximum distance between L and R of MAX_LR_DIST
        - sum of anchor mapq > 0

    """

    if to_file:
        outfile = open(params.OUTPREF + '.resultati.tsv', 'w')
        
    if to_list:
        outlist = []

    header = [
        
        'chromosome',
        'position',
        'strand',
        
        'support_L',
        'support_R',
        'support_total',
        'anchor_mapq_sum',
        
        'TSD',
        
        'position_L',
        'position_R',
        'dist_LR',
    
        'IS_start',
        'IS_end'
        
        ]
    
    if to_file:
        outfile.write('\t'.join(header) + '\n')
    else:
        sys.stdout.write('\t'.join(header) + '\n')
    
    
    for i in range(len(clusters)):
        
        c = clusters[i]
        
        CHROM = params.chromosomes[0]
        
        """ IS info
       
        """
        te_hits = c[2].te_hits['IS6110']
        TE_START = min(te_hits['aligned_positions'])
        TE_END = max(te_hits['aligned_positions'])
        STRAND = '+' if te_hits['strand']['+'] > te_hits['strand']['-'] else '-'
        
    
        """ Support type
         
        """
        SUPPORT_LEFT = len(c[2].breakpoint[0])
        SUPPORT_RIGHT = len(c[2].breakpoint[1])
        SUPPORT_TOTAL = SUPPORT_LEFT + SUPPORT_RIGHT
        
        if SUPPORT_LEFT == 0 or SUPPORT_RIGHT == 0:
            continue
        
        
        """ Position
        Important part: when pooling results from many strains, insertions will count
        as identical when the have the exact same position. By considering only split reads,
        this should work. 
        
        Remove positional outliers due to stray reads.
        """
        
        POS_LEFT = max(remove_outliers(c[2].breakpoint[0])) if SUPPORT_LEFT else 'NA'
        POS_RIGHT = min(remove_outliers(c[2].breakpoint[1])) if SUPPORT_RIGHT else 'NA'
        
        LR_DIST = POS_RIGHT - POS_LEFT
        if LR_DIST > params.MAX_LR_DIST:
            continue

        POS = POS_LEFT  # The right-most base of the left cluster (5' strand) is considered the insertion site
        

            
        """ Anchor read mapping quality
        
        """
        ANCHOR_MAPQS = sum([te_hits['anchor_mapqs'][k] for k in te_hits['anchor_mapqs']])
        if ANCHOR_MAPQS == 0:
            continue
        

        """ Target site duplication:
        If splitreads overlap, extract the overlapping sequence
        """
      
        position_reverse_sr = min(remove_outliers(c[2].breakpoint[1]))
    
        if position_reverse_sr < POS:
            region = [CHROM, position_reverse_sr + 1, POS]
            TSD = consensus_from_bam(region, tmp +'/' + params.OUTPREF + '.bam', [20, 0])
            
        else:
            TSD = 'NA'
        if len(TSD) > 10:
            TSD = 'NA'
     
        outline = [
            CHROM, POS, STRAND, 
            SUPPORT_LEFT, SUPPORT_RIGHT, SUPPORT_TOTAL, ANCHOR_MAPQS,
            TSD,
            POS_LEFT, POS_RIGHT, LR_DIST,
            TE_START, TE_END
                   ]        
        outline = map(str, outline)
        
        if to_file:
            outfile.write('\t'.join(outline) + '\n')
            
        if to_list:
            outlist.append(outline)
            
        else:
            sys.stdout.write('\t'.join(outline) + '\n')
    
    if to_file:    
        outfile.close()
        
    if to_list:
        return outlist
    
