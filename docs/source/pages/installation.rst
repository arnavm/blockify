.. _installation:

Installation
============

blockify runs on Python (>= 3.4) and is installable via ``pip``

.. code-block:: python

   pip install blockify

Development
-----------

To actively develop blockify, clone from GitHub and switch to the
development branch:

.. code-block:: bash

   git clone https://github.com/arnavm/blockify.git
   cd blockify
   git checkout dev

Unit tests are available from the top-level directory:

.. code-block:: python

   python -m unittest tests.test_basic

Two batteries of tests are provided: ``tests.test_basic`` and
``tests.test_advanced``. For routine development, the basic set of tests
should be sufficient. The advanced suite takes much more time and
fetches several large datasets. It is best used when making major
changes to the code.

Disclaimer
----------

Not to be confused with the `similarly-named Spotify plugin <https://github.com/serialoverflow/blockify>`_.