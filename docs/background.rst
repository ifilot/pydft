.. index:: Background

Background
==========

.. contents:: Table of Contents
    :depth: 3

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

In practice, PyDFT starts from atom-centered grids. Each atom receives a radial
Gauss-Chebychev grid and an angular Lebedev grid. The Becke partitioning then
assigns a smooth molecular weight to every atom at every grid point. These
weights add up to one, so an integral over the molecule can be evaluated as a
sum of weighted atomic-grid integrals. This decomposition is implemented by
:class:`pydft.MolecularGrid`, while the per-atom radial/angular grids are
implemented by :class:`pydft.AtomicGrid`.

Hartree potential
-----------------

Electron-electron repulsion is handled by calculating the Hartree potential
by means of solving Poisson's equation. This equation is solved per fuzzy cell,
as detailed in the seminal paper of Becke.

The educational advantage of this approach is that the Coulomb term can be
inspected in stages. PyDFT projects the density in each atomic cell onto real
spherical harmonics, solves the radial Poisson equation for each
spherical-harmonic channel, interpolates the resulting Hartree-potential
coefficients onto the molecular grid, and finally integrates those values with
pairs of basis functions to build the matrix :math:`\mathbf{J}`.

Exchange-correlation functions
------------------------------

:program:`PyDFT` currently supports two exchange-correlation functions:

* LDA: Slater exchange + SVWN5 for the correlation
* PBE: The (standard) Perdew-Burke-Ernzerhof exchange-correlation functional

Both functionals are evaluated on the numerical molecular grid. For LDA, the
energy density depends only on the local electron density :math:`\rho`. For PBE,
the energy density also depends on the gradient invariant
:math:`\sigma = |\nabla\rho|^2`, so PyDFT additionally caches basis-function
gradients and adds the corresponding GGA matrix contribution.
