![](detettore_ad.png)

First running version for testing with different coverages and read lengths.

## Install
Create a conda environment with the required Python packages.

```{bash}
conda env create -f environment.yaml -n detettore6110
```

## Run 
The only required input are reads in fastq format. For the reference genome and the IS target (IS6110) the defaults in the resources folder are used if not stated otherwise. 

Below is the simplest way to run detettore, using IS6110 as a target and the imputed ancestor MTBC0 ([Harrison et al. 2024](https://doi.org/10.1099%2Fmgen.0.001165)) as a reference. Writes to stdout if no output file path (-o) is provided.

```{bash}
detettore6110.py reads_1.fq.gz reads_2.fq.gz
```

## Output
The first line of the output, starting with #CN, is the IS copy number estimated independently of the reference. This is followed by a header and the table containing the IS insertion sites relative to the reference genome. 

At present **only clear split read insertion signatures** are reported. This means that the reference-independent copy number estimate is usually higher than the number of inferred insertion sites, as insertions into complex regions (repeats, SVs, ...) tend to produce more complicated and ambiguous signatures (or none at all). 