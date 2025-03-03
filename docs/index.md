# HLA typing

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
