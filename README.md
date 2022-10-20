# SSBayesHelena

## Overview

Performs the two-step approach for detecting
selection signatures using genomic data from
diploid individuals and biallelic markers developed by
Gianola et al. (2010) and modifications to carry out
inference based on the Laplace approximation and to
incorporate pedigree information using the likelihood
derived by Martínez et al. (2017).

## Features

- Implements relatively simple Bayesian approaches to infer selection signatures via
Wright’s $F_{ST}$ in diploid organisms using data from biallelic molecular markers.
- Two models and two methods to compute de posterior mean of the are available,
resulting in four approximations.
- The basic method is the original one developed by Gianola et al. (2010), the
package also implements a variant using the Laplace approximation instead of
Monte Carlo integration to compute the posterior mean.
- The other method extents the original Bayesian model to infer allele frequencies by
incorporating pedigree information via a modification of the probability mass
function derived by Martínez et al. (2017), for this approach, the posterior mean can
be computed using the Laplace approximation or Monte Carlo integration.
- Allows using arbitrary model hyperparameters and even posing a different set of
hyperparameters for each subpopulation.


## Installation
To install the development version from GitHub repository:
``` r
# install.packages("devtools")
devtools::install_github("edserranoc/SSBayesHelena")
```

## License

The SSBayesHelena package as a whole is licensed under the GPLv3.0. See the 
[LICENSE](LICENSE) file for more details.
