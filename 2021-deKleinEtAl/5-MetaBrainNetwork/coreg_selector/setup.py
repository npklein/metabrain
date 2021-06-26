import setuptools
from setuptools import find_packages
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="coreg-selector",
    version="1.0.0",
    author="Niek de Klein",
    author_email="niekdeklein@gmail.com",
    description="Small package to select data from a brain coreg sqlite database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/npklein/brain_eQTL/tree/master/GeneNetwork/coreg_selector",
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
          'termcolor',
          'scipy',
          'beautifultable']
)
