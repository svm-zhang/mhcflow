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

``` bash
mhcflow --bam NA12046.tumor.so.bam \
  --hla_ref "$PWD/NA12046_class1/finalizer/NA12046.hla.fasta"
  --bed class1.bed \
  --tag abc_v14.uniq \
  --realn_only \
  --outdir "$PWD/NA12046_tumor" \
  --sample NA12046.tumor
```

The `--realn_only` tells `mhcflow` to run only the fishing and realigning
steps using the sample-specific HLA reference obtained in earlier example.

Now you can skip the mapping step (`--skip-map`) in `lohhla`, and directly detect
LOH events.
