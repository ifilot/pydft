# PyDFT

[![build](https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml/badge.svg)](https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml)
[![build](https://github.com/ifilot/pydft/actions/workflows/build_conda.yml/badge.svg)](https://github.com/ifilot/pydft/actions/workflows/build_conda.yml)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pydft/badges/version.svg)](https://anaconda.org/ifilot/pydft)
[![PyPI](https://img.shields.io/pypi/v/pydft?color=green)](https://pypi.org/project/pytessel/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Python based Density Functional Theory code for educational purposes. The
documentation of PyDFT can be found [here](https://ifilot.pages.tue.nl/pydft/).

## Purpose

PyDFT is a pure-Python package for performing localized-orbital DFT calculations
using Gaussian Type Orbitals. PyDFT currently supports LDA and PBE
exchange-correlation functionals. The purpose of PyDFT is mainly to serve as an
educational tool to explain the inner workings of a DFT calculation. This
program is not intended for professional calculations. It is not particularly
fast nor offers a lot of features that more mature open-source of commercial
packages offer. It does offer a unique insight into a working code and a
considerable effort was made in documenting everything.

## Installation

This code depends on a few other packages. To install this code and its
dependencies, run the following one-liner from Anaconda prompt

```bash
conda install -c ifilot pydft pyqint pylebedev pytessel
```

## Usage

### Check the current version

```python
import pydft
print(pydft.__version__)
```

### Performing a simple calculation

```python
from pydft import MoleculeBuilder, DFT

CO = MoleculeBuilder().get_molecule("CO")
dft = DFT(CO, basis='sto3g', verbose=True)
en = dft.scf(1e-4)
print("Total electronic energy: %f Ht" % en)
```

## Community guidelines

* Contributions to PyDFT are always welcome and appreciated. Before doing so,
  please first read the [CONTRIBUTING](CONTRIBUTING.md) guide.
* For reporting issues or problems with the software, you are kindly invited to
  to open a [new issue with the bug label](https://github.com/ifilot/pydft/issues/new?labels=bug).
* If you seek support in using PyDFT, please 
  [open an issue with the question](https://github.com/ifilot/pydft/issues/new?labels=question)
  label.
* If you wish to contact the developers, please send an e-mail to ivo@ivofilot.nl.

## License

Unless otherwise stated, all code in this repository is provided under the GNU
General Public License version 3.