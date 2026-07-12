.. index:: api_reference

API reference
=============

This page collects the public Python objects that are useful when following the
calculation workflow from code. The narrative chapters explain the physical
ideas; the reference below documents the classes and methods that expose those
ideas programmatically.

Public package interface
------------------------

The top-level :mod:`pydft` package re-exports the commonly used classes and
helpers, including :class:`pydft.DFT`, :class:`pydft.MolecularGrid`,
:class:`pydft.AtomicGrid`, ``pydft.Molecule``, and
``pydft.MoleculeBuilder``.

.. automodule:: pydft

SCF driver
----------

.. autoclass:: pydft.DFT
   :members:
   :show-inheritance:

Numerical grids
---------------

.. autoclass:: pydft.MolecularGrid
   :members:
   :show-inheritance:

.. autoclass:: pydft.AtomicGrid
   :members:
   :show-inheritance:

Spherical harmonics
-------------------

.. automodule:: pydft.spherical_harmonics
   :members:

Exchange-correlation functionals
--------------------------------

.. automodule:: pydft.xcfunctionals
   :members:

Internal numerical helpers
--------------------------

These helpers are included because they support educational inspection of the
numerical algorithms. They are not usually needed for ordinary calculations.

.. automodule:: pydft.angulargrid
   :members:

.. automodule:: pydft.stencil
   :members:
