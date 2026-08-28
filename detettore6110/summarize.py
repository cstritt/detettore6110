import argparse
import detettore6110.lib.summaries as summaries


def get_args():
    
    parser = argparse.ArgumentParser(
        description = """
        Combine and summarize detettore6110 results. The following output is generated: 
            - a file combining all anchors
            - a file comining all reference insertions      
            - two presence-absence matrices based on the anchors, one per side
            - a table with meta data for all unique anchors (i.e. rows of the presence-absence matrix)
            - a file with the copy number per sample
        """)
    
    parser.add_argument('-i', dest='path_detettore_results', required=True,
                        help='Path to the folder containing all detettore results.')
    parser.add_argument('-x', dest='samples',
                        help='List of sample names to process, in case the results folder contains more results than desired.')
    parser.add_argument('-o', dest='outpath', 
                        help='Path to the output directory.')
    parser.add_argument('-r', dest='reference', help='Path to reference genome')
    parser.add_argument('-s', dest='seed_len', default=20,
                        help='Length of the anchor part adjacent to the IS that is used to determine unique insertion events')
    parser.add_argument('-mm', dest='max_mismatches', default = 4,
                        help='Maximum number of mismatches allowed in the seed alignment')
    parser.add_argument('-tsd', dest='tsd_len', nargs='+', type=int, default=[0,3,4],
                        help='Alowable length of the target site duplication.')
    
    args = parser.parse_args()
    return args

def main():
    
    # Load arguments
    argumenti = get_args()
    
    include = []
    if argumenti.samples:
        with open(argumenti.samples) as f:
            include = [line.strip() for line in f]
            
    # Define objects
    anchors = summaries.Anchors()
    reference_insertions = summaries.ReferenceInsertions()

    # Get unique anchors and reference insertions
    anchors.process_anchors(20, 4, argumenti, include, min_num_reads=5, min_z_score=-2)
    reference_insertions.create_insertion_dict(argumenti, include)
    
    # Cluster unique anchors to check if there is any redundancy
    anchors.cluster_unique_anchors()

    # Write presence-absence matrices
    anchors.create_matrix('5', argumenti, to_file=True)
    anchors.create_matrix('3', argumenti, to_file=True)
    
    # Write reference insertions
    reference_insertions.write_ref_insertions(argumenti)
    reference_insertions.write_vcf(argumenti, include)

    # Metadata
    anchors.write_metadata(argumenti,reference_insertions)
    anchors.write_ua_sequences(argumenti)
    anchors.write_anchor_map(argumenti)
    
    # Copy numbers
    anchors.estimate_copy_numbers(argumenti.outpath)
    
    # Print stats
    print(f"Total reference insertions: {reference_insertions.stats['total_insertions']}")
    print(f"Unique reference insertion sites: {len(reference_insertions.stats['unique_insertion_sites'])}")
    print(f"Unique anchors 5 prime : {anchors.stats['nr_unique_sites']['5']}")
    print(f"Unique anchors 3 prime : {anchors.stats['nr_unique_sites']['3']}")
    
    
if __name__ == '__main__':
    main()
