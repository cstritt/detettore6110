from Bio import SeqIO


class Read:
    """ 
    Class for reads that reach into the target IS. 
    Storing information for both the anchor part and the
    part mapping against the IS. 
    
    """
    def __init__(self, read_id):
        """
        Initialize a Read object with a read ID.
        
        Parameters
        ----------
        read_id : str
            ID of the read
            
        """
        self.read_id = read_id
    
    def add_coordinates(self, paf_row):
        """
        Add the coordinates of the read to the object, given a paf row.
        
        Parameters
        ----------
        paf_row : list
            A list of strings, where each element is a field from a paf row.
            
        """
        # Anchor part (read part that does not map against IS)
        query_start = int(paf_row[2])
        query_end = int(paf_row[3])
        query_len = int(paf_row[1])
        
        not_mapping = set(range(query_len)) - set(range(query_start, query_end))
        self.anchor_start = min(list(not_mapping))
        self.anchor_end = max(list(not_mapping))
        
            
        # IS part
        self.target_start = int(paf_row[7])
        self.target_end = int(paf_row[8])
        target_len = int(paf_row[6])
        
        # To which side of the IS does the read map? 
        # Assuming that the start and end of the element are the same in the reads 
        # as in the provided sequence           
        if self.target_start == 0:
            self.side = '5'
        elif self.target_end == target_len:
            self.side = '3'
        else:
            self.side = 'NA'  # What happened here?
            
        # Part of the IS covered by the read
        self.is_range = set(range(self.target_start, self.target_end))
        
    
    def add_sequences(self, read):
        
        """
        Add sequences of the anchor and target parts to the read.
        
        Parameters
        ----------
        read : SeqRecord
            The read from which to extract the sequences.
        """
        anchor_seq = read.seq[self.anchor_start:self.anchor_end]        
        target_seq = read.seq[self.target_start:self.target_end]
                
        self.anchor = SeqRecord(
            anchor_seq,
            id=read.id,
            name = read.id,
            description=f'anchor_{self.side}_{self.anchor_start}-{self.anchor_end}'
            )
        
        self.targetpart = SeqRecord(
            target_seq,
            id=read.id,
            name = read.id,
            description=f'target_{self.side}_{self.target_start}-{self.target_end}' 
        )
        

def parse_paf(paf_file, min_anchor_len=20, min_hit_len=20):
    
    """  Traverse IS alignment to extract coordinates of IS and anchor read parts.
 
    Output a dictionary with read IDs, containing info about both the anchor and the IS part. 
     
    Complications:
        - nested insertions
        - close-by insertions
        
    """
    
    read_dict = {}

    with open(paf_file) as f:
        
        for line in f:
            
            fields = line.strip().split('\t')
            query_len = int(fields[1])
            aln_len = int(fields[10])
            
            if (query_len - aln_len) < min_anchor_len:  # anchor not long enough
                continue
            
            if aln_len < min_hit_len:  # IS part not long enought
                continue
            
            read_id = fields[0]
            
            if read_id in read_dict:
                print(f'Warning: {read_id} already in read_dict')
                continue
            
            read = Read(read_id)
            read.add_coordinates(fields)
            read_dict[read_id] = read
    
    return read_dict


def add_seqs_to_read_dict(read_dict, reads, temp_dir, write_fasta=True):
    
    """
    Add sequences from FASTQ files to the read dictionary and optionally write them to FASTA files.

    This function processes a list of FASTQ files, extracting sequences for
    reads present in the provided read dictionary. The sequences are added to
    each read's corresponding entry in the dictionary. Optionally, the function
    can write the sequences to separate FASTA files for each side ('5' and '3').

    Parameters
    ----------
    read_dict : dict
        A dictionary where keys are read IDs and values are Read objects.
    reads : list
        A list of paths to FASTQ files containing the reads.
    temp_dir : str
        The path to the temporary directory where FASTA files will be written.
    write_fasta : bool, optional
        If True, writes the anchor sequences to FASTA files (default is True).

    """

    fasta_out = {
        '5' : [],
        '3' : [], 
        'NA' : []  
        }
    
    for fastq in reads:
        
        print(fastq)
    
        with gzip.open(fastq, "rt") as fastq_handle:
        
            for read in SeqIO.parse(fastq_handle, "fastq"):
                
                if read.id in read_dict:
                    read_dict[read.id].add_sequences(read)
             
                    side = read_dict[read.id].side
                    fasta_out[side].append(read_dict[read.id].anchor)

    if write_fasta:
        for side in ['5', '3']:  
            with open(f'{temp_dir}/anchors.{side}.fasta', 'w') as fasta_handle:
                SeqIO.write(fasta_out[side], fasta_handle, 'fasta')

        
                       