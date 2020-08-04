import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blockify", # Replace with your own username
    version="0.2.0",
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
