# Detect MHC loss of heterozygosity (LOH) event from paired tumor and normal samples

Detecting HLA LOH events has become a popular area of research, particularly
in a clinical setting where patients are treated with immune checkpoint
inhibitor. Homozygous HLA genotypes resulting from LOH may reduce the efficacy
of these treatment.

One popular program to detect HLA LOH, HLALOH, from McGranahan et al. 2017,
bundles the process of realigning reads from the tumor sample within its
workflow. By separating the realignment into an upstream pipeline, the analysis
module can concentrate solely on LOH detection. This separation also avoids
retyping the same normal sample using a potentially different procedure, which
can lead to inconsistencies in sensitivity and introduce technical bias.

In the following section, we provide a walk-through for generating all
necessary results for LOH detection using a paired normal and tumor dataset.
You can use the `--skip-map` option with the `HLALOH` program to bypass its
realignment process and directly perform LOH detection. Alternatively, you can
also use a re-implementation of `HLALOHA` that is available
[here](https://github.com/svm-zhang/lohhla-mod).

## Walk-through

Assuming you have normal and tumor samples with genomic alignments already
generated as `normal.bam` and `tumor.bam`. First, genotype the MHC allele
using the normal sample:

``` bash
mhcflow --bam normal.bam \
  --hla_ref abc_complete.fasta \
  --bed hla.bed \
  --tag abc_uniq.v14 \
  --freq HLA_FREQ.txt \
  --min-ecnt 1 \
  --outdir ./normal
```

The command above produces a sample-level HLA reference and realignment files
in the `finalizer` directory as follows:

``` { .sh }
NA18740_class1
├─ finalizer/
│  ├─ normal.hla.fasta
│  ├─ normal.hla.nix
│  ├─ normal.hla.realn.bam
│  ├─ normal.hla.realn.bam.bai
```

Next, use the sample-level HLA reference to re-align reads from the
paired tumor sample in the `realn-only` mode:

``` bash
mhcflow --bam tumor.bam \
  --hla_ref normal.hla.fasta \
  --bed hla.bed \
  --tag abc_uniq.v14 \
  --freq HLA_FREQ.txt \
  --outdir ./tumor \
  --realn-only
```

???+ note
    The `abc_complete.fasta` reference is replaced with `normal.hla.fasta` to
    ensure that the tumor sample is re-aligned against the reference containing
    the MHC alleles identified from the paired normal sample.

???+ note
    However, the same set of Kmer sequences is used for fishing. This approach
    is beneficial because:<br>
    1. Including Kmer patterns from a larger set of alleles increases
    the fisher's sensitivity.<br>
    2. Using the same set ensures consistent fisher
    sensitivity across different tumor samples.

???+ note
    In the tumor sample command, the population frequency file provide via
    `--freq` is not used--since allele typing is not performed on the tumor
    sample. The `--realn-only` mode ensures to terminate
    `mhcflow` as soon as the `realigner` component completes.

The `mhcflow` command for the tumor sample produces the realignment in the
`realigner` directory:

``` { .sh }
tumor
├─ realigner/
│  ├─ log/
│  ├─ tumor.hla.realn.bam
│  ├─ tumor.hla.realn.bam.bai
```
## Run LOH detection

You can now run [LOH detection](https://github.com/svm-zhang/lohhla-mod)
using the realignments for both normal and tumor samples along with the
sample-level HLA reference:

``` bash
lohhlamod --subject sbj_1 \
  --tbam tumor.hla.realn.bam \
  --nbam normal.hla.realn.bam \
  --hla_ref normal.hla.fasta \
  ...
```


