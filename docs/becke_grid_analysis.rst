Becke grid analysis
===================

.. contents:: Table of Contents
    :depth: 3

All numerical integrations are performed by means of Gauss-Chebychev and Lebedev
quadrature using the Becke grid. 

Atomic fuzzy cells
------------------

It is possible to produce a contour plot of
the fuzzy cells (molecular weights) on a plane using the
:meth:`pydft.MolecularGrid.calculate_weights_at_points` method. An example 
is provided below.

.. literalinclude:: scripts/06-becke-grid.py
	:language: python
	:linenos:
	:emphasize-lines: 22

In the script above, we color every point by the maximum value for each of the
atomic weights. When this maximum value is one, this implies that the grid point
belongs to a single atom. When multiple atoms 'share' a grid point, the maximum
value among the atomic weights will be lower than one.

.. image:: _static/img/user_interface/06-becke-grid.png

.. note::

	Producing such a contour plot is only meaningful for planar molecules such
	as benzene. For more complex molecules such as methane, it is rather
	difficult to make sense of the fuzzy cells upon projection on a plane.

Grid points
-----------

To obtain the set of grid points on which the numerical integration (quadrature),
is performed, we can invoke the :meth:`pydft.MolecularGrid.get_grid_coordinates`
method.

.. literalinclude:: scripts/06-grid-points.py
	:language: python
	:linenos:
	:emphasize-lines: 14

.. image:: _static/img/user_interface/06-grid-points.png

Electron density and its gradient
---------------------------------

The :class:`pydft.DFT` class exposes its internal matrices and a number of
useful functions which we can readily use to interpret its operation. Let us
start by generating a plot of the electron density using the method
:meth:`pydft.DFT.get_density_at_points`.

.. literalinclude:: scripts/01-electron-density.py
    :language: python
    :linenos:
    :emphasize-lines: 25

Running the above script yields the electron density.

.. figure:: _static/img/user_interface/01-electron-density.png

Using largely the same code, we can also readily build the electron density
gradient magnitude using :meth:`pydft.DFT.get_gradient_at_points`.

.. literalinclude:: scripts/02-electron-density-gradient.py
    :language: python
    :linenos:
    :emphasize-lines: 25

.. figure:: _static/img/user_interface/02-electron-density-gradient.png