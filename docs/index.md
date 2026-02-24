# HLA typing

[![PyPI version](https://img.shields.io/pypi/v/mhcflow)](https://pypi.org/project/mhcflow/)
![Python versions](https://img.shields.io/pypi/pyversions/mhcflow)
![License](https://img.shields.io/pypi/l/mhcflow)

## Introduction

`mhcflow` is an end-to-end workflow designed to accurately genotype
MHC class I and II alleles. It streamlines the process of generating
analysis-ready outputs that support both HLA loss of heterozygosity (LOH)
detection and peptide binding prediction.

## Features

`mhcflow` builds upon the well-established Polysolver algorithm,
adding modern enhancements and additional functionality:

- Supports both class I and
  [II](https://github.com/svm-zhang/mhcflow?tab=readme-ov-file#extend-to-class-ii-typing)
  typing

    mhcflow expands on original approach by supporting both class I
  and class II typing while maintaining high accuracy.

- Modern, modular design
    - Re-engineered with a modular architecture for better flexibility and maintainability

    - Streamlined workflows that minimize I/O operations

    - Optimized runtime, ensuring faster analysis without compromising quality

    - Minimal hard-coded logic, facilitating easier customization and integration

- User-friendly integration

    Easily incorporate `mhcflow` into modern whole-genome, whole-exome
    pipelines through its command line interface.

- Analysis-ready

    Generates results can be directly used for HLA LOH detection and
    peptide binding prediction.

## Installation

Starting from `v0.2.0`, `mhcflow` can be installed from PyPI:

```bash
pip install mhcflow
```

If you prefer to use the shell script implementation, please refer to the
instruction on the [Installation](https://svm-zhang.github.io/mhcflow/install) page.


## Citation

Please cite the original
[Polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper.

If you use `mhcflow`, please cite this github repo as well.
