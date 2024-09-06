
# Create versioned direct file and locked dependency file (https://pythonspeed.com/articles/conda-dependency-management/)
environment.yml: versioned main dependencies 
conda-lock.yml: locked subdependencies

```
conda create -n detettore6110
conda activate detettore6110
conda install \
  cd-hit \
  minimap2 \
  pip \
  pysam \
  samtools \
  seqtk \
  pandas

pip install biopython


conda env export --from-history > environment.yml # manually edit to keep main dependencies??

# Update environment
conda env update --file environment.yaml --prune

# Create locked version
conda-lock -f environment.yaml -p linux-64
conda-lock render -p linux-64 # allows using mamba create --file conda-linux-64.lock

# The conda-linux-64.lock file is then used in assemblySC.def to build the container
