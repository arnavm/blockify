.. _tutorial:

Tutorial
========

This tutorial illustrates basic blockify usage. We will analyze some previously published data :cite:`moudgil_self-reporting_2020`. Specifically, we will segment bulk SP1-*piggyBac* calling cards data and call SP1-directed peaks. We will also segment wild-type *piggyBac* calling cards data and call BRD4-directed peaks.

Getting Started
---------------

We need to download and decompress the processed data files, as well as a reference set of TTAA tetramers.

.. code-block:: bash

   > wget -O HCT-116_PBase.ccf.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471636&format=file&file=GSM4471636%5FHCT%2D116%5FPBase%2Eccf%2Etxt%2Egz"
   > gunzip HCT-116_PBase.ccf.gz
   > wget -O HCT-116_SP1-PBase.ccf.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4471637&format=file&file=GSM4471637%5FHCT%2D116%5FSP1%2DPBase%2Eccf%2Etxt%2Egz"
   > gunzip HCT-116_SP1-PBase.ccf.gz
   > wget https://gitlab.com/arnavm/calling_cards/-/raw/master/Ref/TTAA/hg38_TTAA.bed.xz
   > xz -d hg38_TTAA.bed.xz
   > ls
   HCT-116_PBase.ccf
   HCT-116_SP1-PBase.ccf
   hg38_TTAA.bed

Here, HCT-116_SP1-PBase.ccf are the insertions from the SP1-directed experiment, HCT-116_PBase.ccf are the insertions from the wild-type transposase, and hg38_TTAA.bed are the set of potential *piggyBac* insertion sites.

Calling SP1 Peaks
-----------------

We first segment the experiment file as this gives us a candidate set of regions with piecewise-constant densities. These contiguous, bookended intervals are known as blocks.

.. code-block:: bash

   > blockify segment -i HCT-116_SP1-PBase.ccf -o HCT-116_SP1-PBase.blocks
   > wc -l HCT-116_SP1-PBase.blocks
   21375
   > head -n 3 HCT-116_SP1-PBase.blocks
   chr1	54672	758707
   chr1	758707	906326
   chr1	906326	906925

We now use this set of blocks along with the insertions themselves to call peaks. Since the blocks have piecewise-constant density, we model the number of insertions in each block as a Poisson process. For each block, we parameterize the background Poisson process using the control, undirected insertions, scaled for library size. A set of peaks might be called like this:

.. code-block:: bash

   > blockify call -i HCT-116_SP1-PBase.ccf -r HCT-116_SP1-PBase.blocks -bg HCT-116_PBase.ccf -a 0.05 --correction fdr_bh -d 250 -t -o HCT-116_SP1-PBase_peaks.bed
   > wc -l HCT-116_SP1-PBase_peaks.bed
   8356

Here, we specified the input file with ``-i``, the input regions with ``-r``, and the background data with ``-bg``. We used Benjamini-Hochberg correction at a false-discovery rate of 5% (``-a 0.05 --correction fdr_bh``), merging significant blocks within 250 bp of each other (``-d 250``), and tightening the final peaks (``-t``) to improve peak resolution. For a full list of options, see :ref:`call`.

Calling BRD4 Peaks
------------------

To identify BRD4-bound peaks from undirected insertions, we follow a similar set of commands for calling TF-directed peaks. There are two key differences: first, we use hg38_TTAA.bed as the background file, as our null model assumes insertions would be uniformly distributed across the genome; and second, we set the pseudocount to 0 (``-c 0``). When calling TF-directed peaks, both the input and background sets of insertions are random variables. Thus, in any given block, it is possible that there are zero insertions in the background file. To account for undersampling, we added a pseudocount (default: 1). This is not necessary when calling BRD4 peaks because there will always be a TTAA in each block, as *piggyBac* virtually always inserts into this motif. Retaining a pseudocount, while not technically wrong, will decrease sensitivity.

.. code-block:: bash

   > blockify call -i HCT-116_PBase.ccf -r HCT-116_PBase.blocks -bg hg38_TTAA.bed -c 0 -p 1e-30 -d 12500 -o HCT-116_PBase_peaks.bed
   > wc -l HCT-116_PBase_peaks.bed
   1935

Here, we've optimized peak calling to detect BRD4-bound super-enhancers. We've set a very strict threshold, using an unadjusted p-value cutoff of 1e-30 (``-p 1e-30``) and merging significant blocks within 12,500 bp (``-d 12500``).
