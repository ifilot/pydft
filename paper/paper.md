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
Kohn-Sham density functional theory (DFT) with Gaussian-type orbitals. It
implements restricted Kohn-Sham calculations for closed-shell systems and
supports local density approximation (LDA) and generalized-gradient
approximation (GGA) calculations, including the Perdew-Burke-Ernzerhof (PBE)
functional [@pbe:1996]. PyDFT is designed primarily as an educational
implementation: clarity, inspectability, and access to intermediate quantities
are prioritized over production-level efficiency.

PyDFT uses PyQInt [@filot:2025:pyqint] for molecular construction, basis-set
handling, and analytic one-electron integrals. It then constructs the
Kohn-Sham machinery explicitly: atom-centered molecular grids, Becke
fuzzy-cell weights, density and gradient evaluation, Hartree potentials,
exchange-correlation matrices, and the self-consistent field (SCF) procedure.
The code exposes matrices, grid quantities, orbitals, potentials, timing data,
and energy components through a Python interface that can be inspected,
visualized, and modified with standard scientific tools.

# Statement of Need

Density functional theory is one of the central methods of modern electronic
structure theory. Since the Hohenberg-Kohn theorems [@hohenberg:1964] and the
Kohn-Sham formulation [@kohn:1965], DFT has become a standard framework for
computing electronic structures across chemistry, physics, and materials
science. This success also means that many students first encounter DFT through
large research codes. These programs are indispensable, but they are built to
deliver reliable calculations rather than to reveal implementation-level
details. Learners may study the Kohn-Sham equations and exchange-correlation
functionals in a textbook, yet still not see how a molecular DFT calculation is
assembled from numerical grids, density values, potentials, matrix elements,
and convergence controls.

Excellent theoretical treatments exist, including the classic text by Parr and
Yang [@parr:1989]. However, compared with Hartree-Fock theory,
where implementation-level introductions are more common [@szabo], fewer
resources guide readers through the construction of a working molecular DFT
program at the level of code, data structures, grids, and intermediate
matrices. PyDFT addresses this gap. It is not intended to compete with mature
electronic-structure engines, but to provide a compact and modifiable
implementation that shows how the pieces fit together.

A key design goal is that learners can explore implementation choices rather
than treating them as fixed constants. In production software, grid sizes,
Lebedev orders, spherical-harmonic limits, finite-difference stencils, density
normalization, and convergence criteria are often hidden behind defaults or
hardcoded profiles. PyDFT exposes these settings directly, supporting
experiential learning through active experimentation, observation, reflection,
and further testing [@kolb:1984].

# Educational Context

PyDFT has been trialed in both bachelor- and master-level education at
Eindhoven University of Technology within the Chemical Engineering and
Chemistry curriculum. At bachelor level, students use the software as an
accessible calculation environment: they perform DFT calculations, inspect
energies, orbitals, densities, and matrices, and connect numerical results to
chemical structure and bonding. At master level, the emphasis shifts from using
DFT to understanding DFT as an implemented algorithm. Students work more
directly with the source code, molecular grids, Kohn-Sham iteration,
exchange-correlation evaluation, and numerical settings. The same package can
therefore first serve as a transparent calculator and later as an inspectable
model implementation.

PyDFT is used alongside *Elements of Electronic Structure Theory*
[@eoesbook], which presents the quantum-chemical and electronic-structure
theory in systematic form. PyDFT complements this text by making the theory
operational: concepts such as basis functions, one-electron matrices, electron
density, numerical integration, Kohn-Sham orbitals, exchange-correlation
terms, and SCF convergence can be studied in the book and then traced through
executable Python code.

# Features and Implementation

PyDFT implements closed-shell restricted Kohn-Sham DFT in a localized Gaussian
basis. PyQInt supplies molecular and basis-set objects and evaluates the
overlap, kinetic-energy, and nuclear-attraction matrices. PyDFT forms the
one-electron Hamiltonian, orthonormalizes the basis, builds the density matrix,
and iterates the Kohn-Sham matrix to self-consistency with density mixing and
direct inversion of the iterative subspace (DIIS) [@pulay:1980].

The program flow is summarized in Figure 1. The central object is the molecular
grid, based on Becke's multicenter integration scheme [@becke:1988:multicenter].
Each atom receives a radial Gauss-Chebychev grid and an angular Lebedev grid;
Becke fuzzy-cell weights combine the overlapping atom-centered grids into a
molecular quadrature. The Hartree contribution follows the related
Becke-Dickson approach [@becke:1988:poisson]: the density in each fuzzy cell is
projected onto real spherical harmonics, Poisson's equation is solved for the
radial components, and the resulting potential is integrated with Gaussian
basis-function pairs. LDA calculations combine Slater exchange [@slater:1951]
with Vosko-Wilk-Nusair correlation [@vosko:1980], while GGA calculations use
PBE [@pbe:1996] with the density and its gradient.

![High-level overview of the PyDFT workflow. Molecular input, charge, and basis-set information are combined with PyQInt, NumPy, and PyLebedev to construct the numerical objects required by PyDFT. The package then builds atom-centered Becke grids, evaluates densities and potentials, iterates the restricted Kohn-Sham equations to self-consistency, and exposes intermediate matrices, fields, orbitals, energy components, timing data, and tunable settings for educational inspection.](images/pydft-overview.png){#fig:pydft-overview}

The public interface returns both final results and intermediate objects. Users
can inspect matrices, orbital coefficients and energies, density values,
gradients, Becke weights, Hartree and exchange-correlation potentials, energy
components, and timing data. The source code is extensively commented, and the
documentation contains executable examples for visualizing density fields,
orbitals, Becke cells, spherical-harmonic projections, matrix heatmaps, and
energy decompositions. The public interface also exposes the numerical choices
that shape a calculation, including the radial and angular resolution of the
molecular grid, the size of the spherical-harmonic expansion used for the
Hartree potential, the finite-difference treatment of Poisson's equation, the
choice of exchange-correlation approximation, density normalization behavior,
and SCF convergence criteria. Users can therefore study how these
implementation choices affect accuracy, stability, and computational cost.

# Availability

PyDFT is distributed as an open-source Python package and can be installed from
PyPI together with its dependencies. The documentation provides installation
instructions, API documentation, and executable examples for reproducing the
visualizations and numerical experiments described above.

# References
