from setuptools import setup, find_packages

setup(
    name="SHOOT",
    version="1.1.3",
    description="Phylogenetic gene search and ortholog inference",
    author="David Emms",
    author_email="david_emms@hotmail.com",
    url="https://github.com/davidemms/SHOOT",
    packages=["shoot"],
    install_requires=["biopython", "ete3", "numpy", "biopython", "sklearn", "scipy"],
)
