# Behind the scene

`mhcflow` features a modular design and is composed of
four components: `fisher`, `realigner`, `typer`, and `finalizer`. The purpose
of this page is to address the following questions:

- What does each component do?
- Which filter(s) are used (both explicitly and implicitly)?
- What are the limitations of the current implementation?

Refer to the
[diagram](https://svm-zhang.github.io/mhcflow/diagram/#mhcflow-Diagram)
for an overview of `mhcflow` workflow.

## Fisher

The `fisher` component extracts reads sequenced from DNA molecules originating
from MHC class I and II genes. It employs two primary strategies to identify
HLA-derived reads:

- Reads with Kmer sequence matches
- Reads mapped to the HLA regions specified in a BED file

### Kmer-matching reads

`mhcflow` utilizes the Aho-Corasick algorithm (via the
[pyahocorasick](https://pypi.org/project/pyahocorasick/) package) to efficiently
detect reads containing matching Kmer sequences. In this process, only reads
that are mapped to chromosome 6 or are unplaced in the genomic alignments
are analyzed. This targeted (greedy) approach is based on
experiments showing that the vast majority of Kmer-matching reads
originate from chromosome 6, while reads from other chromosomes contribute
negligibly. Due to the polymorphic nature of MHC sequences, `mhcflow`
also examines unplaced reads that have Kmer matches.

???+ note
    Unplaced reads refer to those that do not have a determined placement
    anywhere on a given reference-they are not considered unmapped by definition.

???+ note
    A limitation of this Kmer fishing strategy is that it only captures
    reads with exact matches. Consequently, any sequencing errors reduces
    its efficiency. To compensate for this limitation, 
    `mhcflow` also extracts all reads from a predefined list of
    HLA regions (regardless of the presence of Kmer sequence pattern),
    as described in the next section.
    

### Reads from HLA regions

Extracting all reads from a predefined list of HLA regions enables `mhcflow`
to capture HLA-derived reads that might be missed by the Kmer-fishing process
due to sequencing errors.

???+ note
    In its simplest form, you can specify only the regions to the HLA genes you
    intend to type, which generally produces good results. Alternatively, you
    may specify additional regions to potentially capture more HLA-derived
    reads. However, based on experiments with 1000 genome samples, obtaining more
    fished reads does not necessarily translate to improved typing accuracy.


### Other strategies

`mhcflow` implements the original fishing strategy from the `polysolver`
algorithm, which extracts reads from genomic alignments. There are alternative
approaches to capturing HLA-derived reads. For example, `Optitype`
directly aligns all reads from a sequencing sample against the HLA reference.

- Advantages
    - Captures HLA-derived reads with improved sensitivity. It not only
        detects a greater number of reads but also identify reads that
        are more likely to be HLA-derived.
    - Does not rely on genomic alignment and can be initiated at
        the begining of a workflow after reads trimming.
- Disadvantages
    - HLA-derived reads predominantly originate from
    chromosome 6, as mentioned earlier.
    - The razorS3 aligner used by `Optitype` is not memory-efficient. You
    must find a balance betten the number of threads and
    per-thread memory usage, which can be challenging with varying
    sequencing library sizes.


## Realigner

The realigner component takes the HLA-derived reads and re-aligns them against
the HLA reference using `novoalign`. The realignment process can run in parallel
when `mhcflow` is executed with a `--nproc` value greater than 1.

???+ note
    The following options are used in the `novoalign` command line for
    realignment: `-R 0 -r All -o FullNW`

???+ note
    The realigner component consumes most of the `mhcflow` runtime compared
    to the other components. The number of reads processed by each realignment
    is determined by dividing the total number of HLA-derived reads
    by the value of `--nproc`. <br>
    <br>
    You can profile your workflow to determine the optimal number of reads
    per process based on your assay and computation infrastructure. <br>
    <br>
    For example, in a well-profile assay, the run-to-run variation
    in the number of DNA molecules originating from HLA genes 
    is expected to be minimal. As a result, samples from a study targeting
    the similar coverage are likely to yield a consistent number
    of HLA-derived reads. Due to variation in pull-down efficiency among
    different HLA alleles in targeted assay, it is recommended to determine
    the optimal number of reads by profiling those samples with high counts of
    HLA-derived reads.
    

## Typer

Please refer to the [mhctyper documentaton](https://svm-zhang.github.io/mhctyper/)
for further details on the `typer` component.

???+ note
    Filters applied silently by `mchtyper` can be found
    [here](https://svm-zhang.github.io/mhctyper/#filters-applied-silently)

## Finalizer

The `finalizer` component generates data ready for downstream analyses:

- Sample-level HLA reference: HLA sequences corresponding to the typed
    allele for each sample.
- Sample-level HLA realignment: Realignment of HLA-derived reads against the
    sample-level HLA reference.

???+ note
    In case of presence of homozygous genotypes, the number of reference
    sequences in the sample-level HLA reference may not match the
    number of rows in the typing result table.

???+ note
    The realignment output is coordinate-sorted and indexed but does not
    have duplicates marked.
    
