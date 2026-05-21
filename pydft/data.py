"""
Static data and numerical parameters for atom-centered grids.

This module collects element-specific constants and empirically chosen
numerical parameters used in the construction of atom-centered integration
grids for electronic-structure calculations. The data defined here include:

  • Element symbols and atomic numbers
  • Element-specific radial scale parameters (in Bohr)
  • Recommended numbers of radial shells per atom
  • Conservative angular momentum cutoffs for Lebedev angular grids

The values provided in this file are not meant to represent fundamental
physical constants, but rather well-tested numerical choices that ensure
stable and systematically convergent integration accuracy across a wide
range of chemical elements and environments.

These parameters are intended to be treated as read-only configuration
data. They may be adjusted or extended if higher numerical precision is
required or if alternative grid construction strategies are employed.
"""

ANGSTROM2BOHR = 1.88973

elements = [
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
]

radii = [0.35, 2.00, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50,
         2.25, 1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 2.50,
         2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35,
         1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 2.75,
         2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35,
         1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 3.00,
         2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85,
         1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55,
         1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90,
         1.80, 1.60, 1.90, 1.65, 3.25, 2.80, 2.15, 1.95, 1.80,
         1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
         1.75, 1.75, 1.75, 1.75, 1.55, 1.55]

BSRADII = {
    el: r * ANGSTROM2BOHR
    for el, r in zip(elements, radii)
}

ATCHARGE = {
    el: Z
    for Z, el in enumerate(elements, start=1)
}

# Choice of the number of radial shells per atom (N_shells).
#
# The number of radial shells controls the resolution of the radial
# integration grid and determines how well core, valence, and diffuse
# regions of the electronic density are represented. In practice, the
# required radial resolution increases with the spatial extent and
# complexity of the atomic density.
#
# We therefore choose N_shells based on the principal quantum number of
# the valence shell, grouping atoms by periodic-table rows. Within a
# given row, atoms are assigned the same number of shells, reflecting
# their similar radial structure.
#
# The values used here are conservative and chosen to ensure stable and
# systematically convergent results across chemical environments, rather
# than representing a minimal or variationally optimized choice. This
# strategy is consistent with common practice in real-space and
# atom-centered grid constructions used in electronic-structure methods.
#
# As with the angular grids, these values should be understood as
# numerical parameters controlling integration accuracy, and may be
# increased if higher precision is required.
ATOM_NSHELLS = {
    # First row
    "H":  20,
    "He": 20,

    # Second row (Li–Ne)
    "Li": 25,
    "Be": 25,
    "B":  25,
    "C":  25,
    "N":  25,
    "O":  25,
    "F":  25,
    "Ne": 25,

    # Third row (Na–Ar)
    "Na": 30,
    "Mg": 30,
    "Al": 30,
    "Si": 30,
    "P":  30,
    "S":  30,
    "Cl": 30,
    "Ar": 30,
}

# Choice of angular momentum cutoff l_max for Lebedev angular grids.
#
# Lebedev grids of degree L integrate all spherical harmonics Y_lm exactly
# up to l = L. However, in practical electronic-structure calculations one
# often integrates products of spherical harmonics (e.g. densities, squares
# of orbitals, or nonlinear functionals), which formally require higher
# angular accuracy.
#
# Following the rationale of Becke (J. Chem. Phys. 88, 2547 (1988)), who
# recommended using angular grids whose effective resolution corresponds
# to roughly half the formal angular degree in order to ensure numerical
# stability, we adopt a conservative cutoff
#
#     l_max ≈ L / 2
#
# for Lebedev grids. While Becke’s original discussion was formulated for
# Gaussian-type angular grids, the same reasoning is applied here and
# extrapolated to higher-order Lebedev grids, yielding a robust and
# systematically convergent choice of l_max.
#
# This mapping therefore prioritizes numerical stability over the formal
# maximum angular exactness of the grid.
LMAX_NANGPTS = {
      6:   1,   # L = 3
     14:   2,   # L = 5
     26:   3,   # L = 7
     38:   4,   # L = 9
     50:   5,   # L = 11
     74:   6,   # L = 13
     86:   7,   # L = 15
    110:   8,   # L = 17
    146:   9,   # L = 19
    170:  10,   # L = 21
    194:  11,   # L = 23
    230:  12,   # L = 25
    266:  13,   # L = 27
    302:  14,   # L = 29
    350:  15,   # L = 31
    434:  17,   # L = 35
    590:  23,   # L = 47
    770:  29,   # L = 59
    974:  35,   # L = 71
   1202:  41,   # L = 83
   1454:  47,   # L = 95
   1730:  53,   # L = 107
   2030:  59,   # L = 119
   2354:  65,   # L = 131
   2702:  71,   # L = 143
   3074:  77,   # L = 155
   3470:  83,   # L = 167
   3890:  89,   # L = 179
   4334:  95,   # L = 191
   4802: 101,   # L = 203
   5294: 107,   # L = 215
   5810: 113    # L = 227
}