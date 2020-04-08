rkmh2: Kmer-based read filtering using MinHash
----------------------------------------------

![C/C++ build/test for rkmh2](https://github.com/edawson/rkmh2/workflows/C/C++%20build/test%20for%20kramer/badge.svg)

## Introduction
**rkmh2** is a reimplementation of (some) of the algorithms for MinHash-based read filtering from the
original [rkmh package](https://github.com/edawson/rkmh) and [accompanying paper](https://doi.org/10.1186/s12859-019-2918-y). It focuses on improved stability, lower memory usage, usability and speed.

## Usage
Currently only one command is supported: `rkmh2 filter`

### filter
rkmh filter selects

