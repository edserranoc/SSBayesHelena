# SSBayesHelena

## Overview

 A two-step procedure presented in 
 [Gianola, D., Simianer, H., & Qanbari, S. (2010)](https://www.cambridge.org/core/journals/genetics-research/article/twostep-method-for-detecting-selection-signatures-using-genetic-markers/4447A599402A4EF9862088D3B034B48B) 
 is implemented for analysis of $F_{ST}$ statistics obtained for a battery of loci,
which eventually leads to a clustered structure of values. The first step uses a simple Bayesian model
for drawing samples from posterior distributions of $\theta$ - parameters. This step assigns a weakly informative prior to allelic frequencies and does not make any
assumptions about evolutionary models. The second step regards samples from these posterior
distributions as ‘data’ and fits a sequence of finite mixture models, with the aim of identifying
clusters of $\theta$ - statistics.

## Features

- The first step, drawing the samples without constructing Markov chains. This makes faster the analyses.

- It is implemented a function to compute the Lagrange approximation for the posterior median and variance of $F_{ST}$ parameter. 
## Installation
To install the development version from GitHub repository:
``` r
# install.packages("devtools")
devtools::install_github("edserranoc/SSBayesHelena")
```

## License

The SSBayesHelena package as a whole is licensed under the GPLv3.0. See the 
[LICENSE](LICENSE) file for more details.
