# Getting started

## Quick start

### Input

`mhcflow` requires

- Sorted genomic alignment in BAM format with index
- HLA reference sequence in Fasta and Nix (novoalign index)
- BED file with region of each HLA allele
- HLA kmer tags
- HLA 4-digit supertype frequency table

### Run mhcflow on sample from 1000 genome

Let us type class 1 alleles for `NA12046` sample provided by the 1000
genome project (`NA12046` will be used as example throughout this doc):

``` bash
mhcflow --bam NA12046.so.bam \
  --hla_ref abc_complete.fasta \
  --bed class1.bed \
  --tag abc_v14.uniq \
  --freq HLA_FREQ.txt \
  --outdir "$PWD/NA12046_class1" \
  --sample NA12046
```

## Explain Output

The command above generates HLA typing results in the designated output
directory specified by `--outdir` option. The structure of the output folder
looks like:

``` { .sh }
.
├─ NA12046_class1/
│  ├─ finalizer/
│  ├─ fisher/
│  ├─ realigner/
│  └─ typer/
```

### Fisher output

### Realigner output

### Typer output

`mhcflow` uses `mhctyper` for HLA allele typing. Please refer to
[documentation](https://svm-zhang.github.io/mhctyper/#output-explain) for more details on typing outputs.

### Finalizer output

The `finalizer` folder provides the sample-level HLA reference sequence that
can be directly used as reference for realigning tumor data in a paired tumor
and normal setting. There is also an alignment result in BAM against this
sample-level reference with suffix `hla.realn.ready.bam`. In context of oncology
or immuno-oncology research, the BAM file can be directly ported to program to
detect HLA loss of heterozygosity.

The typing result can be found within the `typer` folder, with suffix
`hlatyping.res.tsv`. The result table should look like the following:

| allele           | gene     | tot_scores   | sample   |
| :--------------- | :------- | :----------- | :------- |
| hla_a_01_01_29   | hla_a    | 2452923.4298 | NA12046  |
| hla_a_02_05_01   | hla_a    | 1766396.924  | NA12046  |
| hla_b_50_01_01   | hla_b    | 1332194.9171 | NA12046  |
| hla_b_57_01_08   | hla_b    | 1134814.4428 | NA12046  |
| hla_c_06_02_01_01| hla_c    | 3020505.4303 | NA12046  |
| hla_c_06_02_01_02| hla_c    | 1519636.0349 | NA12046  |
