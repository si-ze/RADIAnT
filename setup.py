from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="RADIAnT",
    version="0.1.0",
    author="Timothy Warwick, Simonida Zehr",
    author_email="szehr@vrc.uni-frankfurt.de",
    description="A pipeline for processing RNA-DNA interaction data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/si-ze/RADIAnT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "snakemake",
        # other dependencies? e.g., "pandas", "numpy", "biopython"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)