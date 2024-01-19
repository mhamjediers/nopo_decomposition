![StataMin](https://img.shields.io/badge/stata-15-blue) ![issues](https://img.shields.io/github/issues/mhamjediers/nopo_decomposition) ![license](https://img.shields.io/github/license/mhamjediers/nopo_decomposition) ![Stars](https://img.shields.io/github/stars/mhamjediers/nopo_decomposition) ![version](https://img.shields.io/github/v/release/mhamjediers/nopo_decomposition) 

# Stata-Module for Nopo Decompositions, version 0.1'

This Stata `.ado` implements a matching-based decomposition analysis of differences in outcomes between two groups following Ñopo (2008) and some convenient postestimation features. Allowed matching types are (coarsened) exact matching, mutlivariate-distance matching, and propensity-score matching, all of which are internally using [`kmatch`](https://github.com/benjann/kmatch) (Jann, 2017).

## Install

```bash
net install nopo, from("https://raw.githubusercontent.com/mhamjediers/nopo_decomposition/master/")
```

## Documentation

Please consider the help file provided with the package via `help nopo` and the methdological details documented in [this file](https://github.com/mhamjediers/nopo_decomposition/blob/main/te.md).

## References:
Jann, B. (2020). KMATCH: Stata module module for multivariate-distance and propensity-score matching, including entropy balancing, inverse probability weighting, (coarsened) exact matching, and regression adjustment. Statistical Software Components. https://ideas.repec.org//c/boc/bocode/s458346.html

Ñopo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. https://doi.org/10/b6tqwq
