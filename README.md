![](detettore_ad.png)

First running version for testing with different coverages and read lengths.

## Install
Create a conda environment with the required Python packages.
At least Python 3.6 for the new string formatting...

conda create -n detettore6110 cd-hit minimap2 samtools pysam pip
conda activate detettore6110
pip install biopython

conda create -n detettore_dev cd-hit minimap2 samtools pysam pip spyder-kernels=2.5


'''
conda env create -f environment.yaml -n detettore6110
'''

## Run 
The only required input are reads in fastq format. For the reference genome and the IS target the defaults in the resources folder are used if not stated otherwise. Here 

```{bash}
detettore6110.py -f strainX_1.fq.gz strainX_2.fq.gz

```

If no output folder (-o) is provided, this will write to the standard output.


## Output
The first line of the output, starting with #CN, is the IS copy number estimated independently of the reference (see ...). 
This is followed by a header and the table containing the IS insertion sites relative to the reference genome. 