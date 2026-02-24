# MHC class II typing

`mhcflow` implements the `polysolver` algorithm and further extends its
application to genotyping class II MHC alleles. Pre-built
resources for class II typing are provided in the `reference/class2/` folder
of the repository:

``` { .sh }
reference
├─ class2/
│  ├─ class2.ref.fasta
│  ├─ class2.ref.38mer.tag
│  ├─ class2.ref.supertype.freq.tsv
│  └─ class2.chr_prefix_free.bed
```

- `class2.ref.fasta`: The class II HLA reference used with `--hla_ref` option.
- `class2.ref.38mer.tag`: The Kmer sequences generated from the class II
    reference, provided to `--tag` option.
- `class2.ref.supertype.freq.tsv`: The population frequency of class II alleles
    at 4-digit resolution, provided to `--freq` option.
- `class2.chr_prefix_free.bed`: A BED file that specifies the regions of
    class II genes. The coordinates in this file are based on the hg19
    genome build.

???+ note
    If you intent to type a specific set of class II alleles, simply replace
    these file according to the instruction in [Build custom MHC reference]().

???+ note
    `mhcflow` does not provide accompanying `.fai` and `.nix` files along with
    the Fasta file. Please generate them using `samtools faidx` and `novoindex`
    prior to running `mhcflow`.

???+ note
    The class II BED file uses `6` as chromosome 6. If your genome reference
    employs a different naming scheme, please update the file accordingly.
    
## Genotyping class II alleles    

Swapping in the class II resources for the `mhcflow` command demonstrated in
[Getting Started](https://svm-zhang.github.io/mhcflow/getting_start/) guide to
genotype class II alleles:

``` bash
mhcflow --bam [Genomic alignment] \
  --hla_ref class2.ref.fasta \
  --bed class2.chr_prefix_free.bed \
  --tag class2.ref.38mer.tag \
  --freq class2.ref.supertype.freq.tsv \
  --outdir [Path ot output folder]
```

## Benchmark

A preliminary benchmark for class II typing was performed using samples from
1000 Genome project. The results were suprisingly promising and are available
[here](https://github.com/svm-zhang/hla_benchmark).

???+ note
    The class II benchmark was performed using the `v0.1.0` branch of the `mhcflow`
    repository. I am still in the process of testing `mhcflow` on a subset of
    those samples.

## Build custom MHC reference

Coming soon...
