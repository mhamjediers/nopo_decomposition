![version](https://img.shields.io/github/v/release/mhamjediers/nopo_decomposition) ![StataMin](https://img.shields.io/badge/stata-15-blue) ![issues](https://img.shields.io/github/issues/mhamjediers/nopo_decomposition) ![license](https://img.shields.io/github/license/mhamjediers/nopo_decomposition) ![Stars](https://img.shields.io/github/stars/mhamjediers/nopo_decomposition)

# Stata-Module for Nopo Decompositions

This Stata `.ado` implements a matching-based decomposition analysis of differences in outcomes between two groups following Ñopo (2008) and some convenient postestimation features. Allowed matching types are (coarsened) exact matching, multivariate-distance matching, and propensity-score matching, all of which are internally using [`kmatch`](https://github.com/benjann/kmatch) (Jann, 2020).

## Installation

To download and install the module from the Statistical Software Components (SSC):
```bash
ssc install nopo
```

Alternatively, install the current version from GitHub (GitHub files might contain the latest updates that have not been currently published in the SSC):
```bash
net install nopo, from("https://raw.githubusercontent.com/mhamjediers/nopo_decomposition/main/")
```

## Suggested citation:
Sprengholz, Maximilian & Maik Hamjediers, 2024. "NOPO: Stata module to perform Nopo (2008) matching decompositions," [Statistical Software Components S459289](https://ideas.repec.org/c/boc/bocode/s459289.html), Boston College Department of Economics. 

## Documentation

Please consider the help file provided with the package via `help nopo` and the methdological details documented in [this file](https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md).

## Bug Reports

Please either open an issue here or drop us a [mail](mailto:maximilian.sprengholz@hu-berlin.de,maik.hamjediers@hu-berlin.de?subject=[nopo]%20Bug%20Report).

## References
Jann, B. (2020). KMATCH: Stata module module for multivariate-distance and propensity-score matching, including entropy balancing, inverse probability weighting, (coarsened) exact matching, and regression adjustment. Statistical Software Components. https://ideas.repec.org//c/boc/bocode/s458346.html

Ñopo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. https://doi.org/10/b6tqwq
