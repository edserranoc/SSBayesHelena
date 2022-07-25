# SSBayesHelena

## Overview

 A two-step procedure presented in Gianola et al. (2010). is implemented for analysis of FST statistics obtained for a battery of loci,
which eventually leads to a clustered structure of values. The first step uses a simple Bayesian model
for drawing samples from posterior distributions of tetha-parameters, but without constructing Markov
chains. This step assigns a weakly informative prior to allelic frequencies and does not make any
assumptions about evolutionary models. The second step regards samples from these posterior
distributions as ‘data’ and fits a sequence of finite mixture models, with the aim of identifying
clusters of tetha-statistics.

## Features

- The first step uses a simple Bayesian model
for drawing samples from posterior distributions of tetha-parameters, but without constructing Markov
chains.
- The second step regards samples from these posterior
distributions as ‘data’ and fits a sequence of finite mixture models, with the aim of identifying
clusters of tetha-statistics.

## Installation
To install the development version from GitHub repository:
``` r
# install.packages("devtools")
devtools::install_github("edserranoc/SSBayesHelena")
```

## License

The SSBayesHelena package as a whole is licensed under the GPLv3.0. See the 
[LICENSE](LICENSE) file for more details.
