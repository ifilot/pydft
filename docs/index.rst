PyDFT: pure-python density functional theory
============================================

.. image:: https://anaconda.org/ifilot/pydft/badges/version.svg
   :target: https://anaconda.org/ifilot/pydft
.. image:: https://img.shields.io/pypi/v/pydft?color=green
   :target: https://pypi.org/project/pydft/
.. image:: https://gitlab.tue.nl/ifilot/pydft/badges/master/pipeline.svg
   :target: https://gitlab.tue.nl/ifilot/pydft/-/commits/master
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`PyDFT` is a pure-Python package for performing localized-orbital DFT
calculations using Gaussian Type Orbitals. 

.. image:: _static/img/pydft_logo_full_512px.png

:program:`PyDFT` currently supports LDA and PBE
exchange-correlation functionals. The purpose of :program:`PyDFT` is mainly to
serve as an educational tool to explain the inner workings of a DFT
calculation. This program is not intended for professional calculations. It is
not particularly fast nor offers a lot of features that more mature
open-source of commercial packages offer. It does offer a unique insight into
a working code and a considerable effort was made in documenting everything.

.. tip::

   More information on the inner workings of :program:`PyDFT` can be obtained
   from the textbook "Elements of Electronic Structure Theory" (specifically
   chapter 4), which is freely available `via this website <https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/>`_.

:program:`PyDFT` has been developed at the Eindhoven University of Technology,
Netherlands. :program:`PyDFT` and its development are hosted on `github
<https://gitlab.tue.nl/ifilot/pydft>`_.  Bugs and feature
requests are ideally submitted via the `github issue tracker
<https://gitlab.tue.nl/ifilot/pydft/issues>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   background
   user_interface
   api
   community_guidelines

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
