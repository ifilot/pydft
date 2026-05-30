PyDFT: pure-python density functional theory
============================================

.. image:: https://img.shields.io/pypi/v/pydft?color=green
   :target: https://pypi.org/project/pydft/
.. image:: https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml/badge.svg
   :target: https://github.com/ifilot/pydft/actions/workflows/build_pypi.yml
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`PyDFT` is a pure-Python package for performing localized-orbital DFT
calculations using Gaussian Type Orbitals. 

.. image:: _static/img/pydft_logo_full_512px.png
   :width: 360px
   :align: center
   :class: no-frame

:program:`PyDFT` serves as an educational tool that illustrates the inner
workings of a density-functional theory calculation. Currently, :program:`PyDFT`
supports LDA and PBE exchange-correlation functionals. While it is not intended
to replace mature open-source or commercial electronic-structure packages, care
has been taken to achieve reasonable performance within the constraints of a
Python implementation. In addition, a strong emphasis has been placed on code
clarity and comprehensive documentation, providing detailed insight into a fully
working DFT code.

.. tip::

   More information on the inner workings of :program:`PyDFT` can be obtained
   from the textbook "Elements of Electronic Structure Theory" (specifically
   chapter 4), which is freely available `via this website
   <https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/>`_.

:program:`PyDFT` has been developed at the Eindhoven University of Technology,
Netherlands. :program:`PyDFT` and its development are hosted on `GitHub
<https://www.github.com/ifilot/pydft>`_.  Bugs and feature requests are ideally
submitted via the `GitHub issue tracker
<https://www.github.com/ifilot/pydft/issues>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   background
   building_molecules
   electronic_structure_calculations
   orbital_visualization
   becke_grid_analysis
   hartree_potential_analysis
   scalar_field_analysis
   api_reference
   benchmarks
   references
   community_guidelines

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
