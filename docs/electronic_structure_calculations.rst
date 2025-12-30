Calculating electronic energy
=============================

To start, we perform a high-level calculation of the electronic structure
of the carbon-monoxide molecule using the PBE exchange-correlation functional.
To perform this calculation, we first have to construct a
:class:`pydft.DFT` object which requires a molecule as its input. Next, we use
the :meth:`pydft.DFT.scf` routine to start the self-consistent field calculation.

.. literalinclude:: scripts/00-start.py
    :language: python
    :linenos:

Performing this calculation shows that the total electronic energy for this
system corresponds to::

    Total electronic energy:      -111.147096 Ht

Result dictionary
-----------------

The result of an SCF calculation is captured in a Python dictionary object.
This dictionary contains all quantities required to analyze, post-process,
and validate the electronic structure calculation. The data layout is shared
between PyDFT and PyQInt to ensure a consistent interface across methods.

.. list-table:: Description of the data contained in the result dictionary
   :widths: 25 75
   :header-rows: 1

   * - Key
     - Description
   * - :code:`energy`
     - Final converged total electronic energy (Hartree).
   * - :code:`energies`
     - Total electronic energy at each SCF iteration.
   * - :code:`ekin`
     - Electronic kinetic energy,
       :math:`\mathrm{Tr}(\mathbf{T}\mathbf{P})`.
   * - :code:`enuc`
     - Electron-nuclear attraction energy,
       :math:`\mathrm{Tr}(\mathbf{V}\mathbf{P})`.
   * - :code:`erepe`
     - Electron-electron Coulomb (Hartree) energy,
       :math:`\mathrm{Tr}(\mathbf{J}\mathbf{P})`.
   * - :code:`enucrep`
     - Nuclear-nuclear repulsion energy.
   * - :code:`ex`
     - Exchange energy (Hartree-Fock or DFT).
   * - :code:`ec`
     - Correlation energy (DFT only).
   * - :code:`exc`
     - Total exchange-correlation energy.
   * - :code:`orbe`
     - Molecular orbital eigenvalues (orbital energies).
   * - :code:`orbc`
     - Molecular orbital coefficient matrix (AO basis).
   * - :code:`density`
     - One-particle density matrix :math:`\mathbf{P}`.
   * - :code:`fock`
     - Fock matrix :math:`\mathbf{F}`.
   * - :code:`hartree`
     - Hartree (Coulomb) matrix :math:`\mathbf{J}`.
   * - :code:`xc`
     - Exchange-correlation matrix (DFT only, ``None`` for HF).
   * - :code:`overlap`
     - Overlap matrix :math:`\mathbf{S}`.
   * - :code:`kinetic`
     - Kinetic energy matrix :math:`\mathbf{T}`.
   * - :code:`nuclear`
     - Nuclear attraction matrix :math:`\mathbf{V}`.
   * - :code:`hcore`
     - Core Hamiltonian matrix
       :math:`\mathbf{H}_\mathrm{core} = \mathbf{T} + \mathbf{V}`.
   * - :code:`mol`
     - Molecular object defining atoms, geometry, and charge.
   * - :code:`nuclei`
     - Nuclear positions and charges.
   * - :code:`cgfs`
     - List of contracted Gaussian basis functions.
   * - :code:`nelec`
     - Total number of electrons.
   * - :code:`time_stats`
     - Dictionary containing timing information for construction and
       SCF iterations.

Showing the electronic steps
----------------------------

.. literalinclude:: scripts/00-verbose.py
    :language: python
    :linenos:

Executing the script above yields the following output::

	001 | E =  -179.237380 | dE = 0.0000e+00 | 0.0010 s
	002 | E =  -106.678953 | dE = 0.0000e+00 | 0.0810 s
	003 | E =  -117.660767 | dE = 7.2558e+01 | 0.0670 s
	004 | E =  -107.204261 | dE = 1.0982e+01 | 0.0660 s
	005 | E =  -117.396119 | dE = 1.0457e+01 | 0.0610 s
	006 | E =  -117.080704 | dE = 1.0192e+01 | 0.0690 s
	007 | E =  -108.040074 | dE = 3.1541e-01 | 0.0580 s
	008 | E =  -107.453445 | dE = 9.0406e+00 | 0.0620 s
	009 | E =  -110.504979 | dE = 5.8663e-01 | 0.0710 s
	010 | E =  -110.484364 | dE = 3.0515e+00 | 0.0680 s
	011 | E =  -109.605523 | dE = 2.0615e-02 | 0.0680 s
	012 | E =  -111.037173 | dE = 8.7884e-01 | 0.0630 s
	013 | E =  -111.136251 | dE = 1.4317e+00 | 0.0710 s
	014 | E =  -111.147174 | dE = 9.9078e-02 | 0.0730 s
	015 | E =  -111.146826 | dE = 1.0923e-02 | 0.0660 s
	016 | E =  -111.146838 | dE = 3.4816e-04 | 0.0720 s
	017 | E =  -111.146838 | dE = 1.2150e-05 | 0.0660 s
	018 | E =  -111.146838 | dE = 8.0277e-08 | 0.0660 s
	Stopping SCF cycle, convergence reached.

Different exchange-correlation functionals
------------------------------------------

To use a different exchange-correlation functional, we can use the 

.. literalinclude:: scripts/00-xc.py
    :language: python
    :linenos:

which yields the following total electronic energies for the :code:`SVWN5` and
:code:`PBE` exchange-correlation functions::

	SVWN:  -111.14709591483225 Ht
	PBE:  -111.65660426438342 Ht