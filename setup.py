import codecs
import os.path
import setuptools


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blockify", 
    version=get_version("lib/blockify/__init__.py"),
    author="Arnav Moudgil",
    author_email="amoudgil@wustl.edu",
    description="Fast and optimal genome segmentation with Bayesian blocks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/arnavm/blockify",
    package_dir={'': "lib"},
    packages=setuptools.find_packages(where='lib'),
    scripts=["bin/blockify"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable ",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    python_requires='>=3.4',
    install_requires=["numpy", "pandas", "scipy", "statsmodels", "pybedtools"],
    keywords="genomics,segmentation,bayesian",
)
