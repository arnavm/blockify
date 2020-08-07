<!--# blockify-->
[![Build Status](https://travis-ci.org/arnavm/blockify.svg?branch=dev)](https://travis-ci.org/arnavm/blockify
)
[![Documentation Status](https://readthedocs.org/projects/blockify/badge/?version=latest)](https://blockify.readthedocs.io/en/latest/?badge=latest)

Fast, mathematically optimal genome segmentation with Bayesian blocks

## Installation

`pip install blockify`

The earliest recommended version of blockify is 0.2.1.

## Usage

Blockify is available as both a Python library and a command line executable.

To use in Python:
```python
from blockify import annotation
from blockify import segmentation
from blockify import normalization
from blockify import downsampling
```

To use from the command line:
`blockify -h`

For more details, please see the [documentation](https://blockify.rtfd.io).

### Development

To actively develop blockify, clone from GitHub and switch to the
development branch:

```
git clone https://github.com/arnavm/blockify.git
cd blockify
git checkout dev
```

Unit tests are available from the top-level directory:

```
python -m unittest tests.test_basic
```

Two batteries of tests are provided: `tests.test_basic` and
`tests.test_advanced`. For routine development, the basic set of tests
should be sufficient. The advanced suite takes much more time and
fetches several large datasets. It is best used when making major
changes to the code.

## Disclaimer

Not to be confused with the [similarly-named Spotify plugin](https://github.com/serialoverflow/blockify).
