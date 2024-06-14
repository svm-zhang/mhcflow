# polysolver Modern

`polysolvermod` is the original [polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) HLA typing algorithm re-engineered
in modern style. It offers almost all aspects of the original algorithm, adds
new features, runs faster, and is friendly to pipeline integration.

## Features

* Supports both class I and [II](https://github.com/svm-zhang/polysolverMod?tab=readme-ov-file#class-ii-hla-typing) typing with good accuracy
* Generates analysis-ready HLA alignments for HLALOH detection
* Re-engineered in modern style with
  * Modular design
  * Faster runtime with internal optimization
  * Minimum I/O operations
  * Minimum hard-coded code
  * Easy integration to pipeline with packaging and better CLI


## Installation

Please refer to [INSATLL](INSTALL.md) for details.

## Quick Start

`polysolvermod` requires
* Sorted genomic alignment in BAM format with index
* HLA reference sequence in Fasta and Nix (novoalign index)
* BED file with region of each HLA allele
* HLA kmer tags
* HLA 4-digit supertype frequency table

Let us type class 1 alleles for `NA12046` sample provided by the 1000
genome project (`NA12046` will be used as example throughout this doc):

```
polysolvermod --bam NA12046.so.bam \
  --hla_ref abc_complete.fasta \
  --bed class1.bed \
  --tag abc_v14.uniq \
  --freq HLA_FREQ.txt \
  --outdir "$PWD/NA12046_class1 \
  --sample NA12046
```

The command above generates HLA typing results in the designated output
directory specified by `--outdir` option. A quick peek into the result
folder looks like:

```
-- NA12046_class1
   -- finalizer
   -- fisher
   -- realigner
   -- typer
```

## Explain Output

The `finalizer` folder provides the sample-level HLA reference sequence that
can be directly used as reference for realigning tumor data in a paired tumor 
and normal setting. There is also an alignment result in BAM against this
sample-level reference with suffix `hla.realn.ready.bam`. In context of oncology
or immuno-oncology research, the BAM file can be directly ported to program to
detect HLA loss of heterozygosity.

The typing result can be found within the `typer` folder, with suffix
`hlatyping.res.tsv`. The result table should look like the following:

```
allele  gene    tot_scores      sample
hla_a_01_01_29  hla_a   2452923.4298    NA12046
hla_a_02_05_01  hla_a   1766396.924     NA12046
hla_b_50_01_01  hla_b   1332194.9171    NA12046
hla_b_57_01_08  hla_b   1134814.4428    NA12046
hla_c_06_02_01_01       hla_c   3020505.4303    NA12046
hla_c_06_02_01_02       hla_c   1519636.0349    NA12046
```

## Step by Step

`polysolvermod` is re-engineered with a modular design. It generally consists of
4 steps: `fishing`, `realigning`, `typing`, and `realigning` (again). Each module implements basic break-and-continue mechanism, meaning that module finished previously will be automatically skipped. Also it is more friendly to integrate
with pipeline/workflow.

### Fisherman: fishing HLA-relevant reads

The original `polysolver` algorithm fishes HLA-related reads via matching
pre-built kmer (tag) sequence and extracting alignments mapped to regions where
HLA class I allele located. `polysolvermod` follows the same strategy and speeds
it up.

```
fisher --mode faster \
  --tag abc_v14.uniq \
  --bed class1.bed \
  --bam NA12046.so.bam \
  --sample NA12045 \
  --out "$PWD/NA12046_class1/fisher/NA12046.fqs.list.txt
```

The result is plain text file with two lines of fished fastq files 
(paired-end reads).

It is important to note there are other approches to fish HLA-relevant reads.
For instance, `Optitype` aligns trimmed reads against the HLA reference using
`razerS3`. From my experience, direct alignment finds more reads and these
reads tend to align better. However, `razerS3` is not quite memory-efficient,
which in my opinion limits its utility, especially your computing platform
is not unlimited. The approach that the original `polysolver` uses provides
decent fishing result.

### Realigner: realigning fished reads to HLA reference

Next the realigner module aligns the fished reads against the provided
class I HLA reference sequence using `novoalign`, same as the original
`polysolver` program. The difference is the realigner module achieves the
reaignment process in parallel to speed things up a bit.

```
realigner --hla_ref abc_complete.fasta \
  --fqs "$PWD/NA12046_class1/fisher/NA12046.fqs.list.txt \
  --sample NA12046 \
  --out "$PWD/NA12046_class1/realigner/NA12046.hla.realn.so.bam
```

Because the academia version of `novoalign` does not support gzipped fastq
file, this step can take up some disk space depending on the sample
sequencing depth that is HLA-related.



## Class II HLA typing


## Disclaimer

## Citation

Please cite the original [polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper

If you use `polysolvermod`, please cite this github repo as well.