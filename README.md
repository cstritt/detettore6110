![](detettore_ad.png)

First running version for testing with different coverages and read lengths.

## Install
Create a conda environment with the required Python packages.


```{bash}
conda env create -f environment.yaml -n detettore6110
```

## Run 
The only required input are reads in fastq format. For the reference genome and the IS target (IS6110) the defaults in the resources folder are used if not stated otherwise. 

This is thus the simplest way to run detettore, using IS6110 as a target and MTBC0 ([Harrison et al. 2024](https://doi.org/10.1099%2Fmgen.0.001165)) as a reference. If no output file path (-o) is provided, this will write to the standard output.

```{bash}
detettore6110.py reads_1.fq.gz reads_2.fq.gz

```


## Output
The first line of the output, starting with #CN, is the IS copy number estimated independently of the reference. This is followed by a header and the table containing the IS insertion sites relative to the reference genome. 