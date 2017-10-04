# SPHERE
Synthetic PHasE Rate Estimator by single metagenome sequence

[![Build Status](https://travis-ci.com/TaskeHAMANO/SPHERE.svg?token=SzpkyWMHFzySqHiz9qDz&branch=master)](https://travis-ci.com/TaskeHAMANO/SPHERE)

## Requirements

* Python3
* Libraries (automatically installed with pip)
    * Numpy
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

Before using sphere, you have to calculate coverage depth. I recommend to use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Samtools](http://www.htslib.org).

```bash
bowtie2-build
gmap_build -f <Genome_sequence_path> <index_prefix>
```

Next, align with metagenome sequence file with constructed index.

```bash
# for single end sequence
bowtie2 -x <index_prefix> -U <Metagenome_sequence_path> -S <Samfile_path>

# for paired end sequence
bowtie2 -x <index_prefix> -1 <Metagenome_foward_sequence_path> -2 <Metagenome_reverse_sequence_path> -S <Samfile_path>
```

Coverage depth of metagenomic sequence is calculated by samtools. In this procedure, please use -aa option to save zero-coverage depth.

```bash
samtools sort <Samfile_path> -O sam -o <Sorted_samfile_path>
samtools depth <Sorted_samfile_path> -aa > <Coverage_depth_path>
```

### Estimation

Now you can use sphere!! This procedure takes a while, because of MCMC algorithm to estimate parameter.

```bash
sphere_estimate <coverage_depth_path> <Estimated_trend_path>
```

Output file would be formatted below style. In this case, PTR of this coverage depth are estimated as 1.63 with bayesian statistical inteval between 1.425 and 1.854 (95% interval) .

|parameters  | mean                   | se_mean                | sd                     | 2.5%                  | 25%                    | 50%                   | 75%                   | 97.5%               | n_eff| Rhat               | 
|------------|------------------------|------------------------|-----------------------|------------------------|-----------------------|-----------------------|-----------------------|----------------------|------|--------------------|
| y0         | 3.244725098767968      | 0.006946732489113307   | 0.06213346426865469   | 3.0817000534826366     | 3.205858122802362     | 3.249033071864881     | 3.296072440168734     | 3.3377861620809997   | 80.0 | 1.0223190057218234 |
| H          | 0.24395037724188748    | 0.00499450391522601    | 0.0356678922466636    | 0.1772370926811655     | 0.21754851275869802   | 0.24244564021620219   | 0.2698540281492474    | 0.30863477786893434  | 51.0 | 0.9918079063711883 |
| O[0]       | -0.0011292604364484538 | 0.015353624163629178   | 0.1373269893234328    | -0.2489190543412273    | -0.09195473632772605  | 0.00757701632451746   | 0.039094325352131885  | 0.309186272958315    | 80.0 | 0.9878109064054629 |
| O[1]       | -0.6578659976000367    | 0.025342164271510138   | 0.22666720803225243   | -0.9793191816757578    | -0.8306432484045103   | -0.7165325865978199   | -0.50546867319251     | -0.14570802682271736 | 80.0 | 0.9977499893773343 |
| y_raw[0]   | -0.04973399529962034   | 0.10398373177565443    | 0.9300587712178728    | -1.5033498740286158    | -0.9060951825723899   | -0.05411914172013754  | 0.6786712750656294    | 1.5020645285562204   | 80.0 | 1.000081375104273  |
| y_raw[1]   | -0.07461144150213532   | 0.08808293689126122    | 0.7878377381867364    | -1.3919123719143867    | -0.6843507627278305   | -0.06969953035230658  | 0.5466277688425799    | 1.37840622901814     | 80.0 | 0.9900508018100893 |
| ...         | ...                   | ...                    | ...                   | ...                    | ...                   | ...                   | ...                   | ...                  | ...  | ...                |
| y_raw[98]  | -0.10346438387542878   | 0.10125956456017199    | 0.9056930789142926    | -1.546088807723561     | -0.8757685698290372   | -0.11118357874367135  | 0.6153171470693455    | 1.483589846960696    | 80.0 | 1.005083970191178  |
| sigma      | 0.0012902006644375937  | 0.00013023655377049747 | 0.0011648711495445556 | 0.00011198025013377596 | 0.0004905216686313984 | 0.0010456322638122395 | 0.0017293227189403542 | 0.005062722429171181 | 80.0 | 0.9922610203471293 |
| y[0]       | 3.244725098767968      | 0.006946732489113307   | 0.06213346426865469   | 3.0817000534826366     | 3.205858122802362     | 3.249033071864881     | 3.296072440168734     | 3.3377861620809997   | 80.0 | 1.0223190057218234 |
| y[1]       | 3.2468091392413925     | 0.006214325041184842   | 0.05558261290547395   | 3.1128481684614786     | 3.2063177468367687    | 3.2477906740897935    | 3.2905554424430252    | 3.3375675787452277   | 80.0 | 1.0110511376995404 |
| ...         | ...                   | ...                    | ...                   | ...                    | ...                   | ...                   | ...                   | ...                  | ...  | ...                |
| y[99]      | 3.2581941923610374     | 0.007729911737674452   | 0.06913843242205439   | 3.129115053720513      | 3.2121303691394827    | 3.263031313073877     | 3.3004034980324755    | 3.4274739491460435   | 80.0 | 0.9960699243068882 |
| lambda[0]  | 25.818629124172627     | 0.16946753506131154    | 1.5157637135056865    | 21.85011058426157      | 24.73317694070498     | 25.845461848552574    | 27.043541356686266    | 28.18263691879644    | 80.0 | 1.0231601413023468 |
| lambda[1]  | 25.893873240903236     | 0.1508262542917771     | 1.349031029552332     | 22.796712535204318     | 24.8414244264876      | 25.888146958469406    | 26.90625921190773     | 28.200693313933595   | 80.0 | 1.0116704036496502 |
| ...         | ...                   | ...                    | ...                   | ...                    | ...                   | ...                   | ...                   | ...                  | ...  | ...                |
| lambda[99] | 26.17734687735354      | 0.2082167261256452     | 1.8623470146775964    | 23.000017590003985     | 25.057820884156886    | 26.17378629712919     | 27.131612646813274    | 31.12039559041608    | 80.0 | 0.9928753717326023 |
| trend[0]   | 0.004606598858158899   | 0.0008609263577451374  | 0.007203026693149317  | 8.521470759110584e-08  | 0.0005613857162643589 | 0.0014160494018129817 | 0.004936060220532905  | 0.027574102209376568 | 70.0 | 0.9937817183956141 |
| trend[1]   | 0.00583624251289054    | 0.0010044813195353032  | 0.00875686512505755   | 5.162537537107799e-07  | 0.0003454625674718021 | 0.0028125665839377717 | 0.008022150851201396  | 0.0349812261251125   | 76.0 | 0.9913141112566377 |
| ...         | ...                   | ...                    | ...                   | ...                    | ...                   | ...                   | ...                   | ...                  | ...  | ...                |
| trend[99]  | 0.0043215357509366285  | 0.0007978343594683739  | 0.006579110671670231  | 1.322368573454064e-05  | 0.0003434176292956229 | 0.0019790597622193346 | 0.004873902081611361  | 0.024203626700653864 | 68.0 | 0.996006881434185  |
| PTR        | 1.632993540906962      | 0.016186032574963205   | 0.11671914079072053   | 1.425447724623758      | 1.545113376998368     | 1.623998812328479     | 1.7155060855486712    | 1.853910233599092    | 52.0 | 0.9912525097489538 |
| lp__       | 8469.933604676078      | 1.3786272816505014     | 8.609524614444807     | 8455.086062908038      | 8463.807336344065     | 8468.900823716858     | 8476.378901528202     | 8487.448595590806    | 39.0 | 1.0150201166275272 |

The meaning of statistics are followed by Stan. In short,

* mean: expected a posteriori(EAP) of the parameter
* se_mean: standard error
* sd: posterior standard deviation
* 2.5%, 25%, 50%, 75%, 97.5%: posterior percentile
* n_eff: Effective sample size
* Rhat: determination index of convergence

### Utilities

Some of utilities for coverage depth analysis are equipped with sphere.

#### Coverage depth visualization

This scripts are useful if you want to know coverage depth trend before using sphere_estimate.

```bash
sphere_dplot <Coverage_depth_path> <Coverage_depth_plot_path>
```

#### Estimated coverage depth trend visualization

This scripts are useful if you want to know estimated depth trend after estimation.


```bash
sphere_mcplot <Coverage_depth_path> <Estimated_trend_path> <Estimated_trend_plot_path>
```

#### Circular stats

Compuete statistical indexes from coverage depth file.

```bash
sphere_cstats <Coverage_depth_path1> [<Coverage_depth_path2> ...] <stats_path>
```

## License

BSD-3-Clause
