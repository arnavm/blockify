Command: normalize
==================

Calculate normalized rates of events in a .ccf file

Required parameters:

* -i, --input INPUT: input .ccf file
* -o, --output OUTPUT: Output file (bedGraph format)

Optional parameters:

* -r, --regions REGIONS: Regions over which to normalize event counts; should be supplied as a BED file. If not provided, the input file will be segmented using Bayesian blocks. (all options from blockify segment are available)
* -k, --libraryFactor LIBRARYFACTOR: Normalization factor for library size (default: 1000000)
* -l, --lengthFactor LENGTHFACTOR: Normalization factor for the length of regions; used to calculate scaled rates of events per interval (default: None)
