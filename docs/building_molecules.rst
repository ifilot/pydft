.. index:: building_molecules

Building molecules
==================

Molecules can be built either directly using the :code:`Molecule` class or via
a :code:`MoleculeBuilder` convenience routine. 

.. important::

	The :code:`Molecule` and :code:`MoleculeBuilder` classes are not implemented
	in :program:`PyDFT` but are obtained from the :program:`PyQInt` module.

Molecule class
##############

To manually build a molecule, one first need to construct a :code:`Molecule`
object after which one or more atoms can be assigned to the molecule. If no
:code:`unit` is specified, it is assumed that all coordinates are given in
atomic units (i.e. Bohr units).

.. code-block:: python
	:linenos:

	from pyqint import Molecule

	mol = Molecule('co')
	mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
	mol.add_atom('O', 0.0, 0.0, 1.2, unit='angstrom')

MoleculeBuilder class
#####################

Alternatively, a Molecule can be constructed from the :code:`MoleculeBuilder` 
class which uses a name. For more information about the :code:`MoleculeBuilder` 
class, please consult the
`PyQInt documentation on the MoleculeBuilder class <https://pyqint.imc-tue.nl/user_interface.html#using-the-moleculebuilder-class>`_.

.. code-block:: python
	:linenos:

	from pyqint import Molecule

	co = MoleculeBuilder().from_name("CO")