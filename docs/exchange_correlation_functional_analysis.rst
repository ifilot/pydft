Exchange-correlation functional analysis
========================================

.. contents:: Table of Contents
    :depth: 3

The exchange and correlation potentials depend on the electron density and once
the latter has been established, we can readily visualize these potentials.

Visualizing exchange potential
------------------------------

To construct the exchange potential at arbitrary grid points, we can use the
script as shown below which uses the :meth:`pydft.MolecularGrid.get_exchange_potential_at_points`
method.

.. literalinclude:: scripts/08-exchange-potential.py
	:language: python
	:linenos:
	:emphasize-lines: 13

.. image:: _static/img/user_interface/08-exchange-potential.png

Visualizing correlation potential
---------------------------------

In a similar fashion as for the exchange potential, we can use the method
:meth:`pydft.MolecularGrid.get_correlation_potential_at_points` to obtain the
correlation potential field.

.. literalinclude:: scripts/08-correlation-potential.py
	:language: python
	:linenos:
	:emphasize-lines: 13

.. image:: _static/img/user_interface/08-correlation-potential.png