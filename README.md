# mhcflow

<!-- toc -->

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Explain Output](#explain-output)
- [Step by Step](#step-by-step)
    * [Fisherman: fishing HLA-relevant reads](#fisherman-fishing-hla-relevant-reads)
    * [Realigner: realigning fished reads to HLA reference](#realigner-realigning-fished-reads-to-hla-reference)
    * [Typer: typing HLA class I genotype](#typer-typing-hla-class-i-genotype)
- [Realigner: generating analysis-ready HLA typing result](#realigner-generating-analysis-ready-hla-typing-result)
- [Extend to Class II typing](#extend-to-class-ii-typing)
- [Scenario: detecting LOH from paired tumor and normal samples](#scenario-detecting-loh-from-paired-tumor-and-normal-samples)
- [License](#license)
- [Disclaimer](#disclaimer)
- [Citation](#citation)

<!-- tocstop -->

## Introduction

`mhcflow` is the original
[polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) HLA
typing algorithm re-engineered in modern style. It offers almost all aspects
of the original algorithm, adds new features, runs faster, and is friendly to
pipeline integration.

## Features

- Supports both class I and
  [II](https://github.com/svm-zhang/mhcflow?tab=readme-ov-file#extend-to-class-ii-typing)
  typing with good
  [accuracy](https://github.com/svm-zhang/hla_benchmark?tab=readme-ov-file)
- Generates analysis-ready HLA alignments for HLALOH detection
- Re-engineered in modern style with
  - Modular design
  - Faster runtime with internal optimization
  - Minimum I/O operations
  - Minimum hard-coded code
  - Easy integration to pipeline with packaging and better CLI

## Installation

Please refer to [INSATLL](INSTALL.md) for details.

## Quick Start

`mhcflow` requires

- Sorted genomic alignment in BAM format with index
- HLA reference sequence in Fasta and Nix (novoalign index)
- BED file with region of each HLA allele
- HLA kmer tags
- HLA 4-digit supertype frequency table

Let us type class 1 alleles for `NA12046` sample provided by the 1000
genome project (`NA12046` will be used as example throughout this doc):

```bash
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

```text
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

```text
allele  gene    tot_scores      sample
hla_a_01_01_29  hla_a   2452923.4298    NA12046
hla_a_02_05_01  hla_a   1766396.924     NA12046
hla_b_50_01_01  hla_b   1332194.9171    NA12046
hla_b_57_01_08  hla_b   1134814.4428    NA12046
hla_c_06_02_01_01       hla_c   3020505.4303    NA12046
hla_c_06_02_01_02       hla_c   1519636.0349    NA12046
```

## Step by Step

`mhcflow` is re-engineered with a modular design. It generally consists of
4 steps: `fishing`, `realigning`, `typing`, and `realigning` (again). Each
module implements basic break-and-continue mechanism, meaning that module
finished previously will be automatically skipped. Also it is more friendly
to integrate with pipeline/workflow.

The `hlapolysolver.sh` script and `mhcflow` binary (after building the package)
demonstrates each following step, if you are interested.

### Fisherman: fishing HLA-relevant reads

The original `polysolver` algorithm fishes HLA-related reads via matching
pre-built kmer (tag) sequence and extracting alignments mapped to regions where
HLA class I allele located. `mhcflow` follows the same strategy and speeds
it up.

```bash
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

```bash
realigner --hla_ref abc_complete.fasta \
  --fqs "$PWD/NA12046_class1/fisher/NA12046.fqs.list.txt \
  --sample NA12046 \
  --out "$PWD/NA12046_class1/realigner/NA12046.hla.realn.so.bam
```

Because the academia version of `novoalign` does not support gzipped fastq
file, this step can take up some disk space depending on the sample
sequencing depth that is HLA-related.

### Typer: typing HLA class I genotype

The typer module is a complete overhaul of the origial perl scripts
`first_allele_calculations.pl` and `second_allele_calcuations.pl`. The original
`polysolver` algorithm types first and second alleles at a locus in two
separated processes. For each HLA class I allele defined in the HLA reference
sequence, it outputs a plain text file with scores. This generates thousands of
files that creates I/O pressure and make the typing process I/O bound. Also,
it takes about 3-4 script calls to type the first allele and makes it hard to
track when error happens.

The `pyhlatyper` written in this repo tires to improves on all aspects:

1. Typing two alleles with one program call
2. Making typing CPU-bound powered by `polars` and `pysam`
3. Processing alignments to calculate scores in parallel
4. Enabling possibility of typing alleles class II alleles
5. Capturing errors in proper way
6. Free of hard-coded code

```bash
pyhlatyper --freq HLA_FREQ.txt \
  --bam "$PWD/NA12046_class1/realigner/NA12046.hla.realn.so.bam \
  --out "$PWD/NA12046_class1/typer/NA12046.hla_typing.res.tsv
```

One important difference of `pyhlatyper` from the original typing scripts is
that it does not actaully use frequency as prior to calculate posterior scores.
This means the `--race` is always `Unknown`. The choice was made because the race
is usually not a known factor when dealing with real-world data.
I probably will remove the `--race` option from CLI for good in the future.

## Realigner: generating analysis-ready HLA typing result

The original `polysolver` finishes after typing is done. `mhcflow` goes
beyond by providing

1. HLA reference sequence specific to the typed sample
2. Alignment against the sample HLA reference

The reason to have this additional step is to get analysis-ready result.
In oncology and/or immuno-oncology research, one of the questions
people has is to know if there is loss of heterozygosity (LOH) occurring in a tumor
sample. [LOHHLA](https://bitbucket.org/mcgranahanlab/lohhla/src/master/) is the
common go-to algorithm to answer the question. However, `lohhla`, before detecting
any LOH event, goes through realigning both normal and tumor samples, despite typing
has been done for the normal sample. Also realignment, in my opinion, belongs to
pipeline. LOH detection algorithm should be simplified to serve what it is designed
for. To have a clearer picture of what I mean, please refer to [tumor and
normal
scenario](https://github.com/svm-zhang/mhcflow?tab=readme-ov-file#scenario-detecting-loh-from-paired-tumor-and-normal-samples)
below.

The final realignment process splits into two steps.
First to extract and index sample-level HLA reference.

```bash
extractor --hla_ref abc_complete.fasta \
  --sample NA12046 \
  --typeres "$PWD/NA12046_class1/typer/NA12046.hla_typing.res.tsv
  --out "$PWD/NA12046_class1/finalizer/NA12046.hla.fasta

```

Then do the realignment against this new reference.

```bash
realigner \
  --hla_ref "$PWD/NA12046_class1/finalizer/NA12046.hla.fasta
  --fqs "$PWD/NA12046_class1/fisher/NA12046.fqs.list.txt \
  --sample NA12046 \
  --mdup \
  --out "$PWD/NA12046_class1/finalizer/NA12046.hla.realn.ready.bam
```

The `--mdup` option marks PCR duplicates so that when counting coverage during
LOH detection, duplicated reads do not get included. If you want to keep duplicates,
simply not using this option.

## Extend to Class II typing

The original `polysolver` algorithm has been well-known for genotyping Class
I alleles. However, in theory it should also be able to apply to the Class II
case, with certain modification as well as a set of Class II references.

To type the Class II alleles, you only need to swap in the new reference
data, and the CLI is the same as we have shown for the Class I case.

I have done some preliminary benchmark on Class II typing using some samples from
1000 genome project. The result is suprisingly not too shady and can be found
[here](https://github.com/svm-zhang/hla_benchmark).

You can also find all Class II-related reference data within the `reference`
folder in this repo.

## Scenario: detecting LOH from paired tumor and normal samples

Detecting LOH event within the HLA region has been one of the popular subjects
scientists want to look into, especially in a clinical cohort where patients
receive immune checkpoint inhibitor treatment. Homozygous HLA genotypes
decrease the diversity of antigen/neo-antigen the immune system can capture.

`mhcflow` can generate LOH analysis-ready inputs for both tumor and paired
normal samples. Here, I only show how to prepare for the tumor data. You can refer
to upstairs for getting the normal data ready.

Again, let us pretend we have a hypothetical tumor data for sample `NA12046`
(apologize for "giving" this sample a tumor, no harmful damage intented):

```bash
polysolvermod --bam NA12046.tumor.so.bam \
  --hla_ref "$PWD/NA12046_class1/finalizer/NA12046.hla.fasta
  --bed class1.bed \
  --tag abc_v14.uniq \
  --realn_only \
  --outdir "$PWD/NA12046_tumor \
  --sample NA12046.tumor
```

The `--realn_only` tells `mhcflow` to run only the fishing and realigning
steps using the sample-specific HLA reference obtained in earlier example.

Now you can skip the mapping step (`--skip-map`) in `lohhla`, and directly detect
LOH events.

## License

- `mhcflow` respects all LICENSE requirement imposed by the original
  `Polysolver` software, and is licensed under GPL-3.

## Disclaimer

- I, by no means, intend to overtake the origianl idea and implementation
  of `Polysolver` algorithm.
- This repo does not distribute `Polysolver` software, as well as all
  its dependencies such as `novoalign` and `novoindex` under commercial licenses.
- `mhcflow` re-engineered only the HLA typing algorithm. All other
  tools in the `Polysolver` suite was not modified and included in this repo.
- `mhcflow` does not necessarily produce identical result as
  `Polysolver` on typing HLA class I alleles.
- Please interpret result at your own discretion when using
  `mhcflow`.
  [`hla_benchmark`](https://github.com/svm-zhang/hla_benchmark) repo provides
  fundamental assessment of `mhcflow` using 1000 genome data on HLA-A,
  HLA-B, HLA-C, HLA-DQB1, and HLA-DRB1.

## Citation

Please cite the original
[Polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper.

If you use `mhcflow`, please cite this github repo as well.
