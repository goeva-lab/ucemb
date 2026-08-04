# UCEMB: untitled contrastive embedding method benchmark

this repository contains code for running/benchmarking contrastive embedding methods for sc-RNAseq data.

## structure

### [`./py`](./py) / [`./R`](./R): code files for running each method

contains code files for running individual methods.
each file expects to be run:

- within an environment w/ all dependencies present
- with the following arguments provided:
    1. the path to an H5AD file for the background dataset (of shape `(n1, m)`)
    1. the path to an H5AD file for the foreground dataset (of shape `(n2, m)`)
    1. the dimensionality of the output embedding space (an integer, `k`)
    1. the path to specifying where to save an output file, which will be a `.npz`-format archive of a 2-D array of shape `(n1+n2, k)`
    1. a contrastive tuning parameter (this parameter may be ignored if the method in question does not take a contrastive parameter as input)

### [`./data`](./data): data input/output storage

[`gen.py`](./data/gen.py): given an input path for a reference H5AD file, generates background and foreground files of test semi-synthetic perturbation simulation data.
it takes the following arguments:

1. the path to a base H5AD file containing a count matrix and spatial data in the `.obsm["spatial"]` slot
1. the path to an output directory to which output files should be saved
1. a seed to use for pseudo random number generation
1. a comma separated list of `top_n` parameters (e.g. number of genes to perturb)
1. a comma separated list of colon seprated pairs of `x/y` coordinates for circle centers (given as fractions of total width/height)
1. a comma separated list of circle radii
1. a comma separated list of `frac` parameters (e.g. fraction of cells within region of interest to actually perturb)
1. a comma separated list of `lambda` parameters (e.g. mean of poisson distribution used as noise)

### [`./docker`](./docker): definitions for OCI images

contains [`dockerfile`](https://docs.docker.com/reference/dockerfile/) definitions for building containers able to run each method (each matches a file from either).
they all define an entrypoint which runs the corresponding code file.

[`build.py`](./docker/build.py) is a dependency-free python file which can be used to build tar archives of OCI images for each method.
it takes the following arguments:

1. the name/path to a docker-compatible executable to use (e.g. `docker`, `podman`, etc.)
1. the path to an SSH key which can provide access to the [`goeva-lab/spatialPeeler`](https://github.com/goeva-lab/spatialPeeler.git) repo (this is required since it is currently private)
1. the path to the directory containing dockerfile definitions
1. the path to the directory in which docker tar archives should be saved

### [`./cc`](./cc): compute canada job submission

[`submit.py`](./cc/submit.py) is a dependency-free python file which is used to submit jobs to a specific Compute Canada cluster.
it takes the following arguments:

1. the name of the cluster to use to schedule jobs (currently only supports values: `rorqual`, `fir`, and `nibi`)
1. the name of the sponsor account to use for job submission
1. the path to the directory containing outputs of [`gen.py`](./data/gen.py)
1. the path to the directory containing outputs of [`build.py`](./docker/build.py)
1. a comma separated list of `k` parameters (e.g. desired dimensionality of the embedding space)
1. a comma separated list of contrastive parameters
1. optional subset of methods to run (if absent submits jobs for all defined methods)

[`batch.sh`](./cc/batch.sh) is the batch file used by above scheduled SLURM jobs.

## example usage of all scripts

```bash
python ./data/gen.py 'data/'*'.h5ad' 'data/sP_bench/inp' '8' '5,10,20' '0.25:0.75,0.5:0.75,0.75:0.75,0.5:0.5,0.65:0.35' '408.13' '0.3,0.5,0.7,1' '1,3,5' &&
  python ./docker/build.py 'podman' ~/'.ssh/github-'*[^.pub] './docker' './docker/tar' &&
  python ./cc/submit.py 'rorqual' 'def-agoeva' 'data/sP_bench/inp' './docker/tar' '15' '0.2'
```
