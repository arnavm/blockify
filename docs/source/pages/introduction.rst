.. _introduction:

Introduction
============

blockify is a fast and optimal genomic peak caller for one-dimensional data (e.g. BED, qBED, CCF).

The package is built around the Bayesian blocks algorithm :cite:`scargle_studies_2013`, which finds the optimal change points in time series data assuming a Poisson counting process. We also implement a dynamic pruning strategy which achieves linear runtime performance :cite:`killick_optimal_2012`. An interactive notebook demonstrating Bayesian blocks can be found `here <https://observablehq.com/d/d2cafaa7d8c1e018>`_.

While Bayesian blocks was originally developed in the astrophysics community for photon-counting experiments, we find that it has applications in genomics. In particular, we use it to analyze transposon calling cards experiments. Calling cards uses a transposase fused to a transcription factor (TF) to deposit transposons near TF binding sites. Bayesian blocks partitions the genome based on the local density of insertions, which in turn are used to identify peaks and candidate TF binding sites. We have also had success using this algorithm to perfom general-purpose genome segmentation.

.. note::
   Recent papers using calling cards include :cite:`shively_homotypic_2019` and :cite:`liu_quantitative_2020`. For examples of Bayesian blocks in practice, see :cite:`cammack_viral_2020` and :cite:`moudgil_self-reporting_2020`.

blockify is best designed to process qBED files :cite:`moudgil_qbed_2020`, although it will work with BED files.

To get started, please see our :ref:`installation` guide and :ref:`tutorial`.

References
----------

.. bibliography:: ../refs.bib
   :cited:
