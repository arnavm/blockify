Command: call
=============

Call peaks in a .ccf file

Required parameters:

* -i, --input INPUT: input .ccf file
* output: output bed file
* Either:

  * -p, --pValueCutoff PVALUECUTOFF: p-value cutoff (NOTE: This is a straight cutoff and will not take into account multiple hypothesis correction!), OR
  * -a, --alpha ALPHA: alpha for multiple hypothesis correction (must be between 0 and 1) AND
  * --correction CORRECTION: if alpha provided, need to specificity method of multiple hypothesis correction. See `statsmodels.stats.multitest <https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html>`_ for a complete list of choices (default: bonferroni)

* -bg, --background BACKGROUND: Background .ccf file

Optional parameters:

* -r, --regions REGIONS: Regions over which to normalize event counts; should be supplied as a BED file. If not provided, the input file will be segmented using Bayesian blocks. (all options from blockify segment are available)
* --intermediate INTERMEDIATE: Intermediate file to write verbose output (CSV format)
* -d, --distance DISTANCE: Merge features closer than this distance (bp)
* --min MIN: Report peaks larger than this cutoff (bp)
* --max MAX: Report peaks smaller than this cutoff (bp)
* -t, --tight: Shrink peak boundaries to overlap data points
* -c, --pseudocount PSEUDOCOUNT: Pseudocount for background regions (default: 1)
