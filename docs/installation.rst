.. _installation:
.. index:: Installation

Installation
============

PyPi
----

Use Python's built-in virtual environments to isolate your installation. You can
replace `.venv` with an alternative folder.

.. code:: bash

	python3 -m venv .venv
	source .venv/bin/activate
	pip install --upgrade pip
	pip install pydft pyqint pylebedev pytessel mendeleev

Testing
-------

To test that your installation is working, you can run the following snippet
of code

.. code:: python

	import pydft
	print(pydft.__version__)

The version number should be returned.

Simple calculation
------------------

To perform a more involved calculation, one can run the following small
example code:

.. code:: python

	from pydft import MoleculeBuilder,DFT

	co = MoleculeBuilder().get_molecule("CO")
	dft = DFT(co, basis='sto3g')
	en = dft.scf(1e-4)
	print("Total electronic energy: %f Ht" % en)

which should return the following result:

.. code::

	Total electronic energy: -111.147096 Ht