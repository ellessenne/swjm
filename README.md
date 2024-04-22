# Analysis of cohort stepped wedge cluster-randomized trials with non-ignorable dropout via joint modeling

This repository contains additional material for the manuscript titled "Analysis of cohort stepped wedge cluster-randomized trials with non-ignorable dropout via joint modeling", which is currently submitted for publication.
A preprint is available on arXiv.

## Table of Content

The content of this repository is organised as follows:

1. The [`01-example-swjm`](01-example-swjm/) folder contains simulated data and code to analyse a stepped wedge trial with the joint modelling approach that we propose in the manuscript;

1. The `02-example-swjm-individual` folder contains simulated data and code to analyse an [individually-randomised stepped wedge](https://pubmed.ncbi.nlm.nih.gov/30225934/) trial (e.g., without any cluster effect);

1. The [`03-simulation-code`](03-simulation-code/) folder contains code that could be used to reproduce the Monte Carlo simulation studies reported in the paper.
   This code will be uploaded to this repository upon acceptance of the manuscript for publication.

## {simswjm} Package

Note that the simulation code requires the {simswjm} R package, which can be installed from [GitHub](https://github.com/RedDoorAnalytics/simswjm) with the following code:

``` r
# install.packages("devtools")
devtools::install_github("RedDoorAnalytics/simswjm")
```

## License

All code included in this repository is released under the [MIT License](LICENSE.md).
If you find this useful, please cite the manuscript in your work.
