Command: downsample
===================

Downsample a .ccf file in proportion to the value column

Required parameters:

* -i INPUT, --input INPUT: Input .ccf file
* -o, --output OUTPUT: Output file (CCF format)
* -n, --number NUMBER: Number of entries to downsample to (cannot exceed length of input file)

Optional parameters:

* -s SEED, --seed SEED: Random seed
* --naive: Sample every row with equal likelihood
