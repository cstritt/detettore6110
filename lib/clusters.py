import os
import subprocess

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter


class AnchorCluster:

    def __init__(self, cluster_nr, side):
        """
        Initialize an AnchorCluster instance with a cluster number and side. 
        Side refers to whether the reads map to the 5' or 3' side of the IS 
        target.

        Parameters
        ----------
        cluster_nr : int
            The cluster number for this anchor cluster.
        side : str
            The side (either '5' or '3') of the anchor cluster.
        """

        self.cluster_nr = cluster_nr
        self.side = side
        self.cluster_id = f'{side}prime_{cluster_nr}'
        self.reads = []
    
    def add_read(self, read_id, read_dict):
        """
        Add a read to the AnchorCluster instance.
        
        Parameters
        ----------
        read_id : str
            ID of the read to add.
        read_dict : dict
            Dictionary with read IDs as keys and anchor sequences as values.
        """
        anchor_rec = read_dict[read_id].anchor
        self.reads.append(anchor_rec)
    
    def align_anchor_reads(self, temp_dir, args):
        
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
            SeqIO.write(self.reads, fasta_handle, 'fasta')
        
        mafft_cmd = [
            'mafft',
            '--thread', args.cpus,
            '--adjustdirection',
            fasta_path
            ]
        
        subprocess.run(mafft_cmd, check=True, stdout=open(alignment_path, 'w'), stderr=subprocess.DEVNULL)
               
    def get_cluster_consensus(self, temp_dir):
        
        """
        Compute the consensus sequence for the cluster based on alignment.

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
        
        alignment_path = os.path.join(
            temp_dir, f'{self.cluster_id}.aligned.fasta')

        aln = AlignIO.read(open(alignment_path), "fasta")
        aln_smry = AlignInfo.SummaryInfo(aln)
        
        self.nr_reads = len(aln)
        self.aln_len = aln.get_alignment_length()
        
        consensus = ''
        n_sites_with_mismatches = 0
        
        for i in range(self.aln_len):
            col = aln_smry.get_column(i)
            count_missing = col.count('-')
            count_present = self.nr_reads - count_missing
            
            # Check on which side the "tail" of the alignment is
            if i == 0:
                self.depth_start = count_present
            if i == (self.aln_len-1):
                self.epth_end = count_present         
            
            # Get consensus base
            if count_present >= 3:
                counter = Counter(col.replace('-', ''))
                base = counter.most_common(1)[0][0]
                if len(counter) > 1:
                    n_sites_with_mismatches += 1
            else:
                base = '-'
                
            consensus += base
        
        self.consensus = consensus.strip('-').upper()
        self.prop_sites_with_mismatches = round(n_sites_with_mismatches / self.aln_len, 2)
        
    def find_reference_position(self):
        """ Map cluster consensi against reference.  
        
        
        """
        pass
        
        
def cd_hit(fasta_path, output_path):
    
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
        '-d', '0',
        '-c', '0.95',
        '-o', output_path,
        '-sc', '1'
        ]
        
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


def parse_clusters(temp_dir, args):
    
    """
    Parse clusters of anchor sequences and compute their consensus.

    This function processes anchor sequences for both 5' and 3' sides by
    clustering them using `cd-hit-est`. For each cluster, it creates an
    `AnchorCluster` instance, adds reads to it, aligns the reads using MAFFT,
    and computes the consensus sequence.

    Parameters
    ----------
    temp_dir : str
        Path to the temporary directory where intermediate files are stored.
    args : class
        Input arguments containing parameters such as number of CPUs for MAFFT.

    Returns
    -------
    anchor_clusters : dict
        A dictionary containing `AnchorCluster` objects for both 5' and 3' sides,
        keyed by cluster IDs.
    """

    anchor_clusters = {'5':{}, '3':{}}

    for side in anchor_clusters:
        
        clusters = cef.cd_hit(
            f'{temp_dir}/anchors.{side}.fasta', 
            f'{temp_dir}/cd_hit_{side}')
        
        for cluster_id in clusters:
            
            anchor_cluster = cef.AnchorCluster(cluster_id, side)
            
            for read in clusters[cluster_id]:
                anchor_cluster.add_read(read, read_dict)  ### ambiguous read_dict!
                
            anchor_cluster.align_anchor_reads(temp_dir, args)
            anchor_cluster.get_cluster_consensus(temp_dir)
            anchor_clusters[side][cluster_id] = anchor_cluster
        
    return anchor_clusters