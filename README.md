# XP-CLR

Code to compute xp-clr values as per [Chen, Patterson & Reich 2010](https://www.ncbi.nlm.nih.gov/pubmed/20086244).
This implementation was written due to found bugs in the source of the original tool.

## Installation

Clone this git repository into your working directory and `cd`.

```
python setup.py install
```

Also available via conda via the `bioconda` channel:
```
conda install xpclr -c bioconda
```

## Use

This is a python module with a convenience script attached. 
It is designed to run on hdf5 files representing genetic data as generated using [scikit-allel](http://alimanfoo.github.io/2017/06/14/read-vcf.html).
Support is available for VCF and text files in same format as original XPCLR tool, but has not been optimised for this. 
Support to be added for `zarr` format soon.

Interface is under development and may change/break in future.

## Documentation

```
usage: xpclr [-h] --out OUT [--format FORMAT] [--input INPUT]
             [--gdistkey GDISTKEY] [--samplesA SAMPLESA] [--samplesB SAMPLESB]
             [--rrate RRATE] [--map MAP] [--popA POPA] [--popB POPB] --chr
             CHROM [--ld LDCUTOFF] [--phased] [--verbose VERBOSE]
             [--maxsnps MAXSNPS] [--minsnps MINSNPS] [--size SIZE]
             [--start START] [--stop STOP] [--step STEP]

Tool to calculate XP-CLR as per Chen, Patterson, Reich 2010

optional arguments:
  -h, --help            show this help message and exit
  --out OUT, -O OUT     output file
  --format FORMAT, -F FORMAT
                        input expected. One of "vcf" (default), "hdf5", or
                        "txt"
  --input INPUT, -I INPUT
                        input file vcf or hdf5
  --gdistkey GDISTKEY   key for genetic position in variants table of hdf5/VCF
  --samplesA SAMPLESA, -Sa SAMPLESA
                        Samples comprising population A. Comma separated list
                        or path to file with each ID on a line
  --samplesB SAMPLESB, -Sb SAMPLESB
                        Samples comprising population B. Comma separated list
                        or path to file with each ID on a line
  --rrate RRATE, -R RRATE
                        recombination rate per base
  --map MAP             input map file as per XPCLR specs
  --popA POPA           filepath to population A genotypes
  --popB POPB           filepath to population A genotypes
  --chr CHROM, -C CHROM
                        Which contig analysis is based on
  --ld LDCUTOFF, -L LDCUTOFF
                        LD cutoff to apply for weighting
  --phased, -P          whether data is phased for more precise r2 calculation
  --verbose VERBOSE, -V VERBOSE
                        How verbose to be in logging. 10=DEBUG, 20=INFO,
                        30=WARN, 40=ERROR, 50=CRITICAL
  --maxsnps MAXSNPS, -M MAXSNPS
                        max SNPs in a window
  --minsnps MINSNPS, -N MINSNPS
                        min SNPs in a window
  --size SIZE           window size in base pairs
  --start START         start base position for windows
  --stop STOP           stop base position for windows
  --step STEP           step size for sliding windows
```

## File formats

`hdf5` as generated from `vcf` by [scikit-allel](http://alimanfoo.github.io/2017/06/14/read-vcf.html).

Or

`vcf` via [scikit-allel](http://alimanfoo.github.io/2017/06/14/read-vcf.html)

Or

`.geno`
Space delimited file, containing 0/1s with sample haplotypes as columns, and rows as SNPs.

`.map`
Space delimited file. 6 columns: ID, chromosome, Genetic Distance, Position, REF, ALT.

For examples of these files see the [fixture](https://github.com/hardingnj/xpclr/tree/master/fixture) folder used for testing.

## Outputs

`modelL`: The likelihood of the best fitting selection coefficient

`nullL`: The likelihood of the null model (selection coefficient = 0.0).

`sel_coef`: The best fitting selection coefficient

`nSNPs`: Number of SNPs in a window. Suggest to ignore windows where this is small

`nSNPs avail`: When number of SNPs in a window is greater than the maximum specified. If this is consistently much larger than nSNPs, it 
suggests your window size is too large.

`start`, `stop`: start and stop positions of the windows.

`pos_start`, `pos_stop`: actual limits of the used SNPs.

`xpclr`: The log likelihood ratio of the best fitting model vs the null: `2 * (modelL - nullL)`

`xpclr_norm`: Normalized `xpclr` ie `(xpclr - np.nanmean(xpclr))/np.nanstd(xpclr)`
