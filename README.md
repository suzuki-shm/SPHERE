# SPHERE
Synthetic PHasE Rate Estimator by single metagenome sequence

[![Build Status](https://travis-ci.com/TaskeHAMANO/SPHERE.svg?token=SzpkyWMHFzySqHiz9qDz&branch=master)](https://travis-ci.com/TaskeHAMANO/SPHERE)

## Requirements

* Python3
* Libraries (automatically installed with pip)
    * Numpy
    * Scipy
    * Pandas
    * matplotlib
    * PyStan

## How to use

### Install

You can install from PyPi

```bash
pip install sphere
```

or from source

```bash
git clone git@github.com:TaskeHAMANO/SPHERE.git
cd SPHERE
pip install .
```

### Pre-calculation

Before using sphere, you have to calculate coverage depth. Here is a example by [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Samtools](http://www.htslib.org).

#### Building index file for alignment

```bash
bowtie2-build -f <Genome_sequence_path> <index_prefix>
```

#### Aligning sequence file to the index file

Next, align with metagenome sequence file with constructed index.

```bash
# for single end sequence
bowtie2 -x <index_prefix> -U <Metagenome_sequence_path> -S <Samfile_path>

# for paired end sequence
bowtie2 -x <index_prefix> -1 <Metagenome_foward_sequence_path> -2 <Metagenome_reverse_sequence_path> -S <Samfile_path>
```

#### Counting coverage depth

Coverage depth of metagenomic sequence is calculated by samtools. In this procedure, please use -aa option to save zero-coverage depth.

```bash
samtools sort -O sam -o <Sorted_samfile_path> <Samfile_path>
samtools depth -a -a <Sorted_samfile_path> > <Coverage_depth_path>
```

#### (Optional) Select coverage depthes for single genome

If the index file contains multiple sequences, select the depthes with sequence name

```bash
grep <Sequence_name> <Coverage_depth_path> > <Selcted_coverage_depth_path>
```

### (Optional) Compression/Noise reduction by rolling median filter

If ...

* The length of sequence is too big to estimate as you have only small memory,
* The coverage depth is too variable to estimate,

it can be smoothed and compressed by median filter command in sphere.
When you setted strided_length=100 and window_length=100, for each 100 bp, the median of coverage depth is sampled.
Note that when we used too greater number in window_length and average coverage depth is shallow, it became easy to fail in estimation. Thus I recommend to use smaller value than 100 for both parameter.

```bash
sphere_filter -s <strided_length> -w <window_length> <Filtered_coverage_depth_path> <Coverage_depth_path>
```

### Estimation

#### MCMC Algorithm

```bash
sphere_estimate -m vonmises -M sampling <Estimated_result_path> [<Coverage_depth_path1> <Coverage_depth_path2> ...]
```

The meaning of statistics are followed by Stan. In short,

* mean: expected a posteriori(EAP) of the parameter
* se_mean: standard error
* sd: posterior standard deviation
* 2.5%, 25%, 50%, 75%, 97.5%: posterior percentile
* n_eff: Effective sample size
* Rhat: determination index of convergence.
    * If it is greater than 1.1, the MCMC trial doesn't reach to convergence

##### Pros.

* Estimate parameters from complicated model
    * When you try to use mixture models, it'll be better to use MCMC algorithm
* Parameter estimation range is available

##### Cons.

* Takes a long time
* Requires much memories

#### EM Algorithm

```bash
sphere_estimate -m vonmises -M optimizing <Estimated_result_path> <Coverage_depth_path1> [<Coverage_depth_path2> ...]
```

Only maximum likelihood estimates(MLE) were calculated.

##### Pros.

* Fast
* Low memory requirement
* Give us a reasonable estimates for simple model

##### Cons.

* Not good for complicated models
* Only give us a MLE
* High dependence for initialization seed
    * It would be better to estimate multiple seeds. And check the variance of the parameters.

### Utilities

Some of utilities for coverage depth analysis are equipped with sphere.

#### Coverage depth visualization

This scripts are useful if you want to know coverage depth trend before using sphere_estimate.
See `sphere_dplot --help` for more details.

```bash
sphere_dplot <plot_path> <Coverage_depth_path>
```

#### Estimated coverage depth trend visualization

This scripts are useful if you want to know estimated depth trend after estimation.
See `sphere_mcplot --help` for more details.


```bash
sphere_mcplot -m vonmises -M sampling <plot_path> <Coverage_depth_path> <Estimated_result_path>
```

#### Circular stats

Compute directional stats from coverage depth file.
You can calculate mean resultant length of the coverage depth by this tool.
See `sphere_cstats --help` for more details.

```bash
sphere_cstats <stats_path> <Coverage_depth_path1> [<Coverage_depth_path2> ...]
```

#### Pewsey's symmetric test

Compute stats and p-value of Pewsey's symmetric test following [(Pewsey, 2002)](https://www.jstor.org/stable/3316098).
See `sphere_symtest` for more details.

```bash
sphere_symtest <stats_path> <Coverage_depth_path1> [<Coverage_depth_path2> ...]
```

## Q&A

### Fitting parameter by optimizing fails, and it returns message as "LS failed, Hessian reset". What is this?

See [previous issue](https://github.com/facebook/prophet/issues/40) submitted to prophet, which also uses Stan as background estimator.

For this problem, there are a few solutions.

* Use Newton algorithm instead of L-BFGS algorithm. If L-BFGS is rejected, it automatically runs.
    * However, it also likely to achieve local minimum.

* Use different random seed in sphere_estimater comamnd. It can be set by `-ss` option.

* Use MCMC algorithm by `-M sampling` option.

## License

BSD-3-Clause
