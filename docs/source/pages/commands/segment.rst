Command: segment
================

Segment a .ccf file using Bayesian blocks

Required parameters:

* -i, --input INPUT: Input .ccf file
* output: Output file (BED format)

Optional parameters:

* --prior PRIOR: Explicit prior on the number of blocks (*not recommended for general use*)
* --p0 P0: Empirical prior based on a specified false-positive rate; must be between 0 and 1 (default: 0.05)
* --method {OP,PELT}: Segment using the optimal partitioning (OP) or pruned exact linear time (PELT) algorithm (default: PELT)
