.. _background:
.. index:: Background

Background
==========

:program:`PyDFT` is a pure-Python package for performing DFT calculations, 
extending upon the functionality of `PyQInt <https://pyqint.imc-tue.nl/>`_
and leveraging a number of other packages such as 
`PyLebedev <https://github.com/ifilot/pylebedev>`_ and
`PyTessel <https://pytessel.imc-tue.nl/>`_ for quadrature on the unit sphere
and generation of isosurfaces, respectively.

.. tip::

   More information on the inner workings of :program:`PyDFT` can be obtained
   from the textbook "Elements of Electronic Structure Theory" (specifically
   chapter 4), which is freely available `via this website <https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/>`_.

Molecular decomposition
-----------------------

Solving the integrals involved in the electronic structure calculation is handled
by means of numerical integration, also termed quadrature. The quadratures are
solved by decomposing the molecule into so-called "fuzzy" cells as documented
in the work of Becke.

Hartree potential
-----------------

Electron-electron repulsion is handled by calculating the Hartree potential
by means of solving Poisson's equation. This equation is solved per fuzzy cell,
as detailed in the seminal paper of Becke.


Exchange-correlation functions
------------------------------

:program:`PyDFT` currently supports two exchange-correlation functions:

* LDA: Slater exchange + SVWN5 for the correlation
* PBE: The (standard) Perdew-Burke-Ernzerhof exchange-correlation functional