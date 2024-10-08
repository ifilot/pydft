# PyDFT

[![build](https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml/badge.svg)](https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml)
[![build](https://github.com/ifilot/pydft/actions/workflows/build_conda.yml/badge.svg)](https://github.com/ifilot/pydft/actions/workflows/build_conda.yml)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pydft/badges/version.svg)](https://anaconda.org/ifilot/pydft)
[![PyPI](https://img.shields.io/pypi/v/pydft?color=green)](https://pypi.org/project/pydft/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

> [!NOTE]  
> Python based Density Functional Theory code for educational purposes. The
> documentation of PyDFT can be found [here](https://pydft.imc-tue.nl/).

## Purpose

PyDFT is a pure-Python package for performing localized-orbital DFT calculations
using Gaussian Type Orbitals. PyDFT currently supports LDA and PBE
exchange-correlation functionals. The purpose of PyDFT is mainly to serve as an
educational tool to explain the inner workings of a DFT calculation. This
program is not intended for professional calculations. It is not particularly
fast nor offers a lot of features that more mature open-source of commercial
packages offer. It **does** offer a unique insight into a working code and a
considerable effort was made in documenting everything.

> [!TIP]  
> Interested in other **education** quantum chemical codes? Have a look at the
> packages below.
> * [PyQInt](https://github.com/ifilot/pyqint) is a hybrid C++/Python (Cython)
>   code for performing Hartree-Fock calculations. This code contains many
>   relevant features, such a geometry optimization, orbital localization and
>   crystal orbital hamilton population analysis.
> * [HFCXX](https://github.com/ifilot/hfcxx) is a full C++ code for performing
>   Hartree-Fock calculations.
> * [DFTCXX](https://github.com/ifilot/dftcxx) is a full C++ code for performing
>   Density Functional Theory Calculations.

## Features

### Numerical integration using Becke grids

Electron-electron interaction terms (both classical as well as
exchange-correlation) are performed by means of extensive numerical integration
schemes performed over so-called Becke grids. Utilizing these grids, molecular
integrals are decomposed into series of weighted atomic integrals.

![Becke grids](img/becke-grid.png)

### Molecular orbital visualization

PyDFT can be readily used alongside [matplotlib](https://matplotlib.org/stable/)
to make figures of molecular orbitals or density fields.

![Molecular orbitals of CO](img/mo_co.png)

### Extensive output

Internal matrices, e.g. overlap or Hamiltonian matrix, are exposed to the user
and can be readily visualized using specific matrix visualization routines.

![Matrices](img/matrices.png)

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
  to open a [new issue with the bug
  label](https://github.com/ifilot/pydft/issues/new?labels=bug).
* If you seek support in using PyDFT, please [open an issue with the
  question](https://github.com/ifilot/pydft/issues/new?labels=question) label.
* If you wish to contact the developers, please send an e-mail to ivo@ivofilot.nl.

## License

Unless otherwise stated, all code in this repository is provided under the GNU
General Public License version 3.