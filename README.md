# Analysis of cohort stepped wedge cluster-randomized trials with non-ignorable dropout via joint modeling

This repository contains additional material for the manuscript titled "Analysis of cohort stepped wedge cluster-randomized trials with non-ignorable dropout via joint modeling", which is currently submitted for publication.
A preprint is [available on arXiv](https://arxiv.org/abs/2404.14840).

The latest released version of this repository is archived on Zenodo at the following DOI: [10.5281/zenodo.13843867](https://zenodo.org/doi/10.5281/zenodo.13843866)

[![DOI](https://zenodo.org/badge/630868056.svg)](https://zenodo.org/doi/10.5281/zenodo.13843866)

## Table of Content

The content of this repository is organised as follows:

1. The [`01-example-swjm`](01-example-swjm/) folder contains simulated data and code to analyse a stepped wedge trial with the joint modelling approach that we propose in the manuscript;

1. The [`02-example-swjm-individual`](02-example-swjm-individual/) folder contains simulated data and code to analyse an [individually-randomised stepped wedge](https://pubmed.ncbi.nlm.nih.gov/30225934/) trial (e.g., without any cluster effect);

1. The [`03-simulation-code`](03-simulation-code/) folder contains code that could be used to reproduce the Monte Carlo simulation studies reported in the paper.
   This code will be uploaded to this repository upon acceptance of the manuscript for publication.

1. Finally, the [`04-convergence-issues`](04-convergence-issues/) folder contains code and steps that could be used to improve convergence of the joint model.

## {simswjm} Package

Note that the simulation code requires the {simswjm} R package, which can be installed from [GitHub](https://github.com/RedDoorAnalytics/simswjm) with the following code:

``` r
# install.packages("remotes")
remotes::install_github("RedDoorAnalytics/simswjm")
```

## License

All code included in this repository is released under the [MIT License](LICENSE.md).
If you find this useful, please cite the manuscript and/or this repository in your work.
