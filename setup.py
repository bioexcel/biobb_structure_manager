"""
Setup module for standard Pypi installation.
"""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biobb_structure_manager",
    version="3.0.2",
    author="Biobb developers",
    author_email="josep.gelpi@bsc.es",
    description="BioBB_structure_manager is a library to efficiently load and process biomolecular 3D structures.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Bioinformatics Workflows BioExcel Compatibility",
    url="https://github.com/bioexcel/biobb_structure_manager",
    project_urls={
        "Documentation": "http://biobb_structure_manager.readthedocs.io/en/latest/",
        "Bioexcel": "https://bioexcel.eu/"
    },
    packages=setuptools.find_packages(exclude=['docs', 'test']),
    install_requires=['biopython==1.76'],
    python_requires='==3.7.*',
    classifiers=(
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ),
)
