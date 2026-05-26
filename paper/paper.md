---
title: 'PyDFT: A Teaching-Oriented Density Functional Theory Implementation in Python'
tags:
  - Quantum chemistry
  - Electronic structure theory
  - Density Functional Theory
  - Gaussian basis functions
  - Molecular integrals
  - SCF method
authors:
  - name: I.A.W. Filot
    orcid: 0000-0003-1403-8379
    corresponding: true
    affiliation: 1
affiliations:
 - name: Inorganic Materials and Catalysis, Department of Chemical Engineering and Chemistry, Eindhoven University of Technology
   index: 1
date: 8 March 2026
bibliography: paper.bib
---

# Summary

PyDFT is a pure-Python package for learning and exploring molecular
Kohn-Sham density functional theory (DFT) using Gaussian-type orbitals. The
package implements restricted Kohn-Sham calculations for closed-shell molecular
systems and currently supports local density approximation (LDA) and
generalized-gradient approximation (GGA) exchange-correlation functionals,
including the Perdew-Burke-Ernzerhof (PBE) functional [@pbe:1996]. PyDFT is
designed primarily as an educational implementation: numerical efficiency and
feature breadth are secondary to clarity, inspectability, and direct access to
the computational steps that are normally hidden inside production electronic
structure packages.

The program uses PyQInt [@filot:2025:pyqint] for molecular construction,
Gaussian basis functions, and the analytic evaluation of overlap, kinetic
energy, and nuclear attraction integrals. PyDFT then builds the remaining
Kohn-Sham machinery explicitly, including atom-centered molecular quadrature
grids, Becke fuzzy-cell weight functions, electron-density and gradient
evaluation, Hartree potentials obtained from a spherical-harmonic solution of
Poisson's equation, exchange-correlation matrices, and the self-consistent field
(SCF) procedure. The package exposes intermediate matrices, grid quantities,
orbital coefficients, orbital energies, potentials, timing information, and
energy components through a Python interface that can be inspected and modified
with standard scientific tools such as NumPy and Matplotlib.

PyDFT is intended for students, instructors, and researchers who want to see how
the abstract equations of molecular DFT become a working program. The source
code is extensively commented, and the documentation contains executable
examples that visualize density fields, Kohn-Sham orbitals, Becke grids,
spherical-harmonic projections, exchange and correlation potentials, and the
effect of changing numerical parameters. In this way, PyDFT functions as both a
calculator for small illustrative systems and a readable computational text on
the implementation of molecular DFT.

# Statement of Need

Density functional theory is one of the central methods of modern electronic
structure theory. Since the Hohenberg-Kohn theorems [@hohenberg:1964] and the
Kohn-Sham formulation [@kohn:1965], DFT has become a standard framework for
computing electronic structures and molecular properties across chemistry,
physics, and materials science. In practical molecular and materials modeling,
generalized-gradient approximations such as PBE occupy a particularly important
role because they combine broad applicability with comparatively modest
computational cost [@pbe:1996].

This success has also made DFT a method that many students encounter first as
users of established software packages. Such packages are indispensable for
research, but they are generally designed to deliver reliable calculations at
scale rather than to reveal implementation-level details. As a result, learners
may see the Kohn-Sham equations, exchange-correlation functionals, and total
energy expressions in a textbook, but still have only a partial view of the
numerical choices that make a molecular DFT calculation possible: how a
multicenter molecular integral is decomposed into atom-centered quadratures, how
the electron density is evaluated on a grid, how the Hartree potential is
constructed, how exchange-correlation potentials become matrix elements, or how
grid parameters affect accuracy and cost.

There are excellent theoretical treatments of DFT, with the classic book of
Parr and Yang [@parr:1989] providing one of the canonical introductions to the
density-functional formulation of atoms and molecules. However, compared with
Hartree-Fock theory, for which implementation-level introductions are more
common in computational chemistry education [@szabo], there are fewer resources
that guide a reader through the construction of a working molecular DFT program
at the level of code, data structures, numerical grids, and intermediate
matrices. PyDFT aims to fill this gap. It does not seek to compete with mature
open-source or commercial DFT engines. Instead, it provides a compact,
transparent, and modifiable implementation that shows how the pieces fit
together.

An important design goal is that learners can explore the consequences of
implementation choices rather than treating them as fixed constants. In many
production programs, choices such as radial grid size, angular quadrature order,
spherical-harmonic expansion limit, density normalization, finite-difference
stencil size, and convergence criteria are hidden behind default settings or
hardcoded profiles. PyDFT exposes these parameters directly. This supports
experiential learning, where understanding develops through cycles of active
experimentation, observation, conceptual interpretation, and further testing
[@kolb:1984]. For example, a student can increase the number of radial shells,
change the Lebedev angular grid, vary the maximum angular momentum in the
spherical-harmonic expansion, or adjust the finite-difference stencil used in
the Poisson solver, and then observe the resulting changes in total energy,
timing, density normalization, and numerical stability.

# Features and Implementation

PyDFT implements closed-shell restricted Kohn-Sham DFT in a localized Gaussian
basis. At the start of a calculation, molecular geometry and Gaussian basis
functions are supplied through PyQInt, which also evaluates the overlap matrix
`S`, kinetic-energy matrix `T`, and nuclear-attraction matrix `V`. PyDFT
constructs the one-electron Hamiltonian `H = T + V`, orthonormalizes the basis
through the eigensystem of `S`, iteratively builds the density matrix `P`, and
forms the Kohn-Sham matrix from one-electron,
electron-electron, and exchange-correlation contributions. The SCF procedure
uses density mixing and direct inversion of the iterative subspace (DIIS)
[@pulay:1980] to accelerate convergence.

The central numerical structure in PyDFT is the molecular grid. Its design is
closely based on the multicenter integration scheme introduced by Becke
[@becke:1988:multicenter]. Each atom receives an atom-centered spherical grid.
The radial coordinate is discretized using a Gauss-Chebychev mapping, while the
angular coordinate is sampled with Lebedev quadrature. To convert the collection
of overlapping atom-centered grids into a molecular quadrature, PyDFT evaluates
Becke weight functions. These smooth fuzzy-cell functions form a partition of
unity over molecular space and allow a molecular integral to be written as a
sum of single-center atomic-like contributions. The implementation exposes both
the generated grid points and the Becke weights, making it possible to visualize
atomic cells, inspect regions where cells overlap, and examine how molecular
integration is assembled from local quadratures.

The Hartree contribution is calculated through a second piece of Becke-inspired
machinery. Following Becke and Dickson [@becke:1988:poisson], PyDFT projects the
electron density in each atom-centered fuzzy cell onto real spherical harmonics.
For every radial shell, the density is represented by coefficients
`rho_klm`, where `k` indexes the radial grid and `l,m` index the spherical
harmonics. Poisson's equation is then solved separately for each
spherical-harmonic component using a finite-difference representation of the
radial differential operator. The resulting Hartree-potential coefficients are
interpolated with cubic splines and evaluated back on the molecular grid. This
procedure gives the Coulomb matrix by numerical integration of the Hartree
potential with pairs of Gaussian basis functions. PyDFT exposes the
spherical-harmonic coefficients, grid potentials, and local atomic-grid objects,
allowing users to inspect not only the final Coulomb matrix but also the
intermediate representation used to construct it.

Exchange and correlation are evaluated numerically on the same molecular grid.
For LDA calculations, PyDFT combines Slater exchange [@slater:1951] with the
Vosko-Wilk-Nusair correlation functional [@vosko:1980]. For GGA calculations,
PyDFT implements the PBE exchange-correlation functional [@pbe:1996], using both
the electron density and its gradient. The code exposes exchange and correlation
energies separately, constructs the corresponding matrix contributions, and
provides helper methods for evaluating exchange and correlation potentials at
arbitrary points in space. This makes the exchange-correlation part of the
calculation visible as an object of study rather than merely as an input keyword
to an SCF driver.

After convergence, the `DFT.get_data()` method returns a structured dictionary
containing the matrices and quantities used in the calculation, including
`S`, `T`, `V`, `J`, `XC`, `F`, `P`, the coefficient matrix
`C`, orbital energies, total and per-iteration energies,
exchange-correlation components, nuclear repulsion, and timing data. The
`DFT.get_molgrid_copy()` method gives access to the molecular-grid object,
through which users can examine grid coordinates, Becke weights, electron
density, density gradients, spherical-harmonic projections, Hartree potentials,
orbital amplitudes, and exchange-correlation potentials. The documentation uses
these methods to build visual examples of density and gradient fields,
Kohn-Sham orbital contour plots, isosurfaces, Becke cells, grid-point
distributions, matrix heatmaps, and energy decompositions.

Several numerical options are intentionally part of the public interface. The
number of radial shells (`nshells`), the number of angular Lebedev points
(`nangpts`), the spherical-harmonic expansion limit (`lmax`), the
finite-difference stencil size (`fdpts`), the exchange-correlation functional,
the density-normalization behavior, convergence tolerance, and optional
parallel construction of spherical harmonics can all be adjusted. The
documentation includes examples that vary these values and report their effect
on total energy and computational time. These controls are not meant primarily
as performance-tuning conveniences; they are part of the pedagogical surface of
the package, allowing users to connect numerical analysis, chemical accuracy,
and computational cost in concrete calculations.

# Availability

PyDFT is distributed as an open-source Python package and can be installed from
PyPI together with its dependencies. The documentation provides installation
instructions, API documentation, and executable examples for reproducing the
visualizations and numerical experiments described above. The project is hosted
on GitHub, where users can report issues, inspect the source code, and
contribute improvements.

# References
