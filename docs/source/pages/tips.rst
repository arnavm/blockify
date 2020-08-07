.. _tips:

Tips and Tricks
===============

Here are some helpful tips and tricks to get the most out of blockify.

Resolution
----------

Transposition data can be sparse, particularly if the transposase is constricted to specific motifs (e.g. *piggyBac*). Sparse data lead to broader peaks, which can be harder to interpret. Here are some strategies to increase resolution:

* Omit or decrease the ``-d/--distance`` flag to minimize merging of significant blocks
* Specify ``-t/--tight``, which will pull in peak boundaries so they overlap qBED entries
* For the sharpest intervals, use the ``-s/--summit`` flag to return each peak's maximum
* Finally, increasing the value of ``--p0`` (default: 0.05) can lead to more peaks being called, at the risk of returning more false positives.

Miscellaneous
-------------

* Although the :ref:`tutorial` demonstrated first generating a list of blocks for input into ``blockify call``, this step is not strictly necessary. If a regions file is not supplied, blockify will generate one behind the scenes using the default settings in ``blockify segment``. However, this can result in considerable memory usage. Pre-computing the blocks file is one way to minimize memory consumption and improve performance.
* Similarly, the regions over which to run ``blockify call`` need not be Bayesian blocks. The program can operate on any set of intervals provided in BED format. This flexibility can be useful if there are a set of features that are biologically meaningful to your analysis. For example, this could be a file of promoter regions or accessible loci where a TF might be bound.
* Peaks are output in BED6 format with a generic annotation, like ``peak_1743``. The program does not re-calculate p-values on these *post hoc* by default. If you want to further calculate significance or density of these peaks, simply re-run ``blockify call`` with the ``--intermediate`` flag set and supply the peaks file to ``-r/--regions``. Then inspect the intermediate file for these details. Picking up with the BRD4 example from the :ref:`tutorial`:

.. code-block::

   > blockify call -i HCT-116_PBase.ccf -r HCT-116_PBase_peaks.bed -bg hg38_TTAA.bed -c 0 -p 1e-30 -d 12500 --intermediate HCT-116_PBase_peaks_annotated.csv > /dev/null
   > head -n 2 HCT-116_PBase_peaks_annotated.csv
   ,chrom,start,end,name,score,strand,Input,Background,Normed_bg,Net_density,pValue,negLog10pValue,rejected
   0,chr1,7298597,7304456,peak_0,1,.,130.0,32,2.533130940448789,0.021755738020063357,3.74243334103237e-169,168.42684592642368,True
