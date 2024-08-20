#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:57:12 2024

@author: cristobal
"""



def ISmap(queries, targets, outpref, min_aln_len, k, w):
    """
    

    Parameters
    ----------
    queries : TYPE
        DESCRIPTION.
    targets : TYPE
        DESCRIPTION.
    outpref : TYPE
        DESCRIPTION.
    min_aln_len : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.

    Returns
    -------
    minimapd : TYPE
        DESCRIPTION.

    """

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


class readCluster:
    """
    Container for read clusters supporting non-reference TE insertion. Here
    mapping (against IS) results are combined with the anchor reads.
    
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



class mise_en_place:
    """ 
    Object for easy access to program settings and input files. 
    Pass args from argparse when initiating object.

    """

    def __init__(self, args):

        # Files
        self.fastq = [os.path.abspath(x) for x in args.fastq]
        self.ref = os.path.abspath(args.ref)
        self.target = os.path.abspath(args.target)

        self.outpref = args.outpref
        self.cpus = args.cpus
        
        # Thresholds
        self.mapq = args.mapq
        self.min_len = args.min_len
        self.max_tsd_len = args.max_tsd_len

        # IS target library
        self.telib = {
            seq_record.id:
                seq_record.seq for seq_record in SeqIO.parse(self.target, "fasta")}
            
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
            
        # Reference
        self.chromosomes = [seq_record.id for seq_record in SeqIO.parse(self.ref, 'fasta')]

        
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