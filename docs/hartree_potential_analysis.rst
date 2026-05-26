.. index:: hartree_potential_analysis

Hartree potential analysis
==========================

.. contents:: Table of Contents
    :depth: 3

The Hartree potential is calculated from the electron density distribution
using the Poisson equation :cite:p:`becke:1988:poisson`.

PyDFT solves this problem in an atom-centered representation. The molecular
density is first split into Becke fuzzy cells. Within each cell, the angular
dependence is expanded in real spherical harmonics, leaving a set of radial
equations that can be solved numerically. The resulting coefficients are later
interpolated back to the full molecular grid to assemble the Coulomb matrix.

Projection onto spherical harmonics
-----------------------------------

The Poisson equation is solved on the **atomic grid** by first projecting 
the electron density onto the spherical harmonics such that the following
equation holds.

.. math::

	\rho(r_{k}, \Omega) = \sum_{lm} \rho_{klm} Y_{lm}

where :math:`\rho(r_{k}, \Omega)` is the electron density at radius
:math:`r_{k}`, :math:`\rho_{klm}` the linear expansion coefficient at position
:math:`r_{k}` for the spherical harmonic :math:`lm` and :math:`Y_{lm}` a
spherical harmonic with quantum numbers :math:`l` and :math:`m`.

We can readily visualize and interpret this projection using the
:meth:`pydft.MolecularGrid.get_rho_lm_atoms` as shown by the script below.

.. literalinclude:: scripts/07-spherical-harmonics-projection.py
	:language: python
	:linenos:
	:emphasize-lines: 13

.. figure:: _static/img/user_interface/07-spherical-harmonics-projection.png
   :width: 78%

   Spherical-harmonic expansion coefficients of the atom-centered density.

Radial Poisson solve
--------------------

After the density coefficients :math:`\rho_{klm}` have been constructed, PyDFT
solves a finite-difference form of the radial Poisson equation for each
:math:`(l,m)` channel. The solution is stored as Hartree-potential expansion
coefficients :math:`U_{klm}` on the radial grid of each atom.

During an SCF calculation the density changes at every iteration, so these
coefficients also change. The geometry of the interpolation problem, however,
does not change: the molecular grid points remain at the same distances from
each atomic center. PyDFT therefore precomputes cubic interpolation stencils once
and reuses them when the updated :math:`U_{klm}` values are needed on the
molecular grid.

Coulomb matrix assembly
-----------------------

Once the Hartree potential is known on the molecular grid, the Coulomb matrix is
obtained by numerical quadrature over products of basis functions:

.. math::

   J_{ij} = \int \chi_i(\mathbf{r})\,U(\mathbf{r})\,\chi_j(\mathbf{r})\,
            d\mathbf{r}.

In the code, this final step is performed by
:meth:`pydft.MolecularGrid.calculate_coulombic_matrix`.
