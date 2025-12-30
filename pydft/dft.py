# -*- coding: utf-8 -*-

from .moleculargrid import MolecularGrid
from pyqint import PyQInt, Molecule, CGF
import numpy as np
import time
from copy import deepcopy
from collections.abc import Mapping
from .data import ATOM_NSHELLS, LMAX_NANGPTS

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 3
SUBSPACE_START = 4

class DFT():
    def __init__(self, 
                 mol:Molecule, 
                 basis:str|list[CGF] = 'sto3g', 
                 functional:str = 'svwn5',
                 nshells: Mapping[str, int] | None = None,
                 nangpts: Mapping[str, int] | None = None,
                 lmax: Mapping[str, int] | None = None,
                 fdpts:int=7,
                 normalize:bool = True):
        """
        Construct a density-functional theory (DFT) calculation object.

        Parameters
        ----------
        mol : Molecule
            Molecular system defining atoms, geometry, and total charge.
        basis : str or list[CGF], optional
            Atomic orbital basis set. This can be specified either as the name
            of a built-in basis set or as an explicit list of contracted
            Gaussian functions. The default is ``'sto3g'``.
        functional : str, optional
            Exchange-correlation functional. Valid options are
            :code:`svwn5` and :code:`pbe`. The default is ``'svwn5'``.
        nshells : Mapping[str, int], optional
            Number of radial integration shells per atomic species. If not
            provided, default values are used based on the atomic row.
        nangpts : Mapping[str, int], optional
            Number of angular integration points per atomic species. The values
            must correspond to supported Lebedev grid sizes. If not provided,
            default values are used.
        lmax : Mapping[str, int], optional
            Maximum angular momentum quantum number per atomic species used in
            the spherical harmonic expansion. If not provided, :code:`lmax` is
            determined automatically from the number of angular points.
        fdpts : int, optional
            Number of grid points used in the finite-difference scheme. The
            default is 7.
        normalize : bool, optional
            Whether to perform intermediate normalization of the electron
            density during the self-consistent field procedure. The default
            is ``True``.
        """
        self.__mol = mol
        self.__calculate_grid_settings(nshells, nangpts, lmax)
        self.__integrator = PyQInt()
        self.__basis = basis
        self.__time_stats = {}
        self.__itermax = 100
        self.__fdpts = fdpts
        self.__functional = functional
        self.__normalize = normalize
        
        # keep track of time
        self.calctimes = {
            'density_hartree': [],
            'calculate_J': [],
            'calculate_XC': [],
            'ulm_interpolation': [],
            'build_hartree_field': [],
            'build_repulsion_matrix': [],
        }

    def get_data(self) -> dict:
        """
        Return results of the SCF calculation.

        This method is only valid after a successful SCF run. It returns a
        dictionary containing all relevant physical quantities, matrices,
        energies, and timing information. The returned data layout is shared
        between PyQInt and PyDFT calculations to ensure consistent
        post-processing.

        Returns
        -------
        dict
            Dictionary containing SCF results and metadata.

        Notes
        -----
        The returned dictionary contains the following entries:

        **System information**
            * ``mol`` : Molecular object defining geometry and atoms
            * ``nuclei`` : Nuclear positions and charges
            * ``cgfs`` : Contracted Gaussian basis functions

        **Energies**
            * ``energy`` : Final total electronic energy
            * ``energies`` : SCF energy history
            * ``ekin`` : Electronic kinetic energy
            * ``enuc`` : Electron-nuclear attraction energy
            * ``enucrep`` : Nuclear-nuclear repulsion energy
            * ``ex`` : Exchange energy
            * ``ec`` : Correlation energy
            * ``exc`` : Exchange-correlation energy

        **Orbital quantities**
            * ``orbc`` : Molecular orbital coefficient matrix
            * ``orbe`` : Molecular orbital eigenvalues

        **Matrices and operators**
            * ``overlap`` : Overlap matrix
            * ``kinetic`` : Kinetic energy matrix
            * ``nuclear`` : Nuclear attraction matrix
            * ``hcore`` : Core Hamiltonian matrix (T + V)
            * ``density`` : Density matrix
            * ``fock`` : Fock matrix
            * ``hartree`` : Hartree (Coulomb) matrix
            * ``xc`` : Exchange-correlation matrix (DFT only, ``None`` for HF)

        **Timing information**
            * ``time_stats`` : Dictionary with timing breakdowns
                - ``construct`` : Grid and setup construction times
                - ``scf`` : SCF iteration timings
        """
        data = {
            # system
            "mol": self.__mol,
            "nuclei": self.__nuclei,
            "cgfs": self.__cgfs,
            "nelec": self.__nelec,

            # core results
            "energy": self.__energies[-1],
            "energies": self.__energies,

            # orbital information
            "orbc": self.__C,        # MO coefficients
            "orbe": self.__e,        # MO eigenvalues

            # density & operators
            "density": self.__P,
            "fock": self.__F,
            "overlap": self.__S,
            "kinetic": self.__T,
            "nuclear": self.__V,
            "hcore": self.__T + self.__V,

            # electron interaction terms
            "hartree": self.__J,
            "xc": self.__XC,

            # energies (explicit)
            "ex": self.__Ex,
            "ec": self.__Ec,
            "exc": self.__Exc,
            "enucrep": self.__enuc,
            "ekin": np.einsum("ij,ji", self.__T, self.__P),
            "enuc": np.einsum("ij,ji", self.__V, self.__P),
            "erepe": np.einsum("ij,ji", self.__J, self.__P),

            # timing
            "time_stats": {
                "construct": self.__molgrid.construct_times,
                "scf": self.calctimes,
            }
        }
        
        return data
    
    def get_molgrid_copy(self) -> MolecularGrid:
        """
        Get a copy of the underlying molgrid object

        Returns:
            MolecularGrid: copy of the :class:`MolecularGrid` object
        """
        return deepcopy(self.__molgrid)
    
    def get_density_at_points(self, spoints:np.ndarray) -> np.ndarray:
        """
        Get the electron density at provided points

        Parameters
        ----------
        spoints : ndarray
            Set of sampling points (:math:`X \\times 3` array)

        Returns
        -------
        ndarray
            Electron density scalar field (:math:`N \\times 1` array)
        """
        if len(spoints.shape) != 2:
            raise Exception('Grid points need to be supplied as a Nx3 array.')

        return self.__molgrid.get_density_at_points(spoints, self.__P)
    
    def get_gradient_at_points(self, spoints:np.ndarray) -> np.ndarray:
        """
        Get the electron density gradient at provided points

        Parameters
        ----------
        spoints : ndarray
            Set of sampling points (:math:`X \\times 3` array)

        Returns
        -------
        ndarray
            Electron density gradient vector field (:math:`N \\times 3` array)
        """
        if len(spoints.shape) != 2:
            raise Exception('Grid points need to be supplied as a Nx3 array.')

        return self.__molgrid.get_gradient_at_points(spoints, self.__P)
    
    def scf(self, tol:float=1e-5, verbose:bool=False) -> dict:
        """
        Perform the self-consistent field procedure

        Parameters
        ----------
        tol : float, optional
            electronic convergence criterion, by default 1e-5

        Returns
        -------
        float
            total electronic energy (in Hartrees)
        """
        # construct stagnant matrices
        self.__setup()
        
        # create empty P matrix as initial guess
        self.__P = np.zeros_like(self.__S)

        # start SCF iterative procedure
        nitfin = 0
        ediff = 0
        for niter in range(0, self.__itermax):
            start = time.perf_counter()
            energy = self.__iterate(niter, 
                                    giis=True if nitfin == 0 else False,
                                    mix=0.9)
            self.__energies.append(energy)
            stop = time.perf_counter()
            itertime = stop - start
            self.__time_stats['iterations'].append(itertime)
            
            if verbose:
                print('%03i | E = %12.6f | dE = %5.4e | %0.4f s' % (niter+1, energy, ediff, itertime))
            
            if niter > 0:
                ediff = np.abs(energy - self.__energies[-2])
            if niter > 2:
                if ediff < tol:
                    # terminate giis self-convergence and continue with mixing
                    nitfin += 1
                    if nitfin < 2:
                        continue
                    
                    # terminate self-convergence cycle
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                        
                        # update density matrix from last found coefficient matrix
                        self.__P = self.__calculate_P()
                    break

        return self.get_data()
    
    def print_time_statistics(self):
        """
        Print a summary of time statistics
        """
        print('-- Construction times --')
        print('Atomic grids:                                %.4f s' % self.__molgrid.construct_times['atomic_grids'])
        print('Projection cart. coord. on solid angles:     %.4f s' % self.__molgrid.construct_times['cartesian_solid_angle_projection'])
        print('Fuzzy cell decomposition:                    %.4f s' % self.__molgrid.construct_times['fuzzy_cell_decomposition'])
        print('Spherical harmonics:                         %.4f s' % self.__molgrid.construct_times['spherical_harmonics'])
        print('Nuclear distance and potential:              %.4f s' % self.__molgrid.construct_times['nuclear_distance_and_potential'])
        print('Basis set amplitudes:                        %.4f s' % self.__molgrid.construct_times['basis_function_amplitudes'])
        print('Spline construction:                         %.4f s' % self.__molgrid.construct_times['spline_construction'])
        print()
        print('-- Calculation times --')
        print('Classical e-e repulsion matrix (J):          %.4f s' % np.average(self.calctimes['calculate_J']))
        print('  - Ulm interpolation:                       %.4f s' % np.average(self.calctimes['ulm_interpolation']))
        print('  - Build Hartree-field:                     %.4f s' % np.average(self.calctimes['build_hartree_field']))
        print('  - Build repulsion matrix:                  %.4f s' % np.average(self.calctimes['build_repulsion_matrix']))
        print('Building edens (rho) and Hartree pot (U):    %.4f s' % np.average(self.calctimes['density_hartree']))
        print('Exchange-correlation matrices (XC):          %.4f s' % np.average(self.calctimes['calculate_XC']))
    
    def get_construction_times(self) -> dict:
        """
        Return construct times dictionary from MolecularGrid to get insights
        into the time to construct particular properties

        Returns
        -------
        dict
            dictionary of construction times for the MolecularGrid objects
        """
        return self.__molgrid.construct_times
    
    def __calculate_grid_settings(self, nshells, nangpts, lmax):
        """
        Build a Mapping where the grid settings for each atom are set, unless
        they are already provided by the user.

        The lmax mapping can be inferred from nangpts, but the user may override
        this.
        """
        # collect atom types
        attypes = np.unique([a[0] for a in self.__mol])

        # build dictionary of number of radial points (Gauss-Chebychev grid)
        # per atom type
        if nshells is None:
            self.__nshells = {}
            for a in attypes:
                self.__nshells[a] = ATOM_NSHELLS[a]
        else:
            self.__nshells = nshells

        # build dictionary of number of angular points (for Lebedev grid) per
        # atom type
        if nangpts is None:
            self.__nangpts = {}
            for a in attypes:
                self.__nangpts[a] = 50 if a == 'H' else 110
        else:
            self.__nangpts = nangpts

        # set lmax values based on number of angular points
        if lmax is None:
            self.__lmax = {}
            for k,v in self.__nangpts.items():
                self.__lmax[k] = LMAX_NANGPTS[v]
        else:
            self.__lmax = lmax

    def __iterate(self, niter, giis=True, mix=0.9):
        """
        Perform single-step iteration
        """
        # calculate J and XC matrices based on the current electron
        # density estimate as captured in the density matrix P
        if niter > SUBSPACE_START and giis:
            try:
                diis_coeff = self.__calculate_diis_coefficients(self.__evs_diis)
                self.__F = self.__extrapolate_fock_from_diis_coefficients(self.__fmats_diis, diis_coeff)
                self.__Fprime = self.__X.transpose().dot(self.__F).dot(self.__X)
                self.__e, self.__Cprime = np.linalg.eigh(self.__Fprime)
                self.__C = self.__X.dot(self.__Cprime)
                self.__P = self.__calculate_P()
            except np.linalg.LinAlgError: # set giis to False if not working
                giis = False
        
        # calculate J and XC matrices based on the current electron
        # density estimate as captured in the density matrix P
        if np.any(self.__P):
            st = time.perf_counter()
            self.__molgrid.build_density(self.__P, normalize=self.__normalize)
            self.calctimes['density_hartree'].append(time.perf_counter() - st)
            
            st = time.perf_counter()
            self.__J = self.__calculate_J()
            self.calctimes['calculate_J'].append(time.perf_counter() - st)
            
            st = time.perf_counter()
            self.__XC, self.__Exc = self.__calculate_XC()
            self.calctimes['calculate_XC'].append(time.perf_counter() - st)
            

        # calculate Fock matrix
        self.__F = self.__H + self.__J + self.__XC

        # perform unitary transformation on Fock matrix
        self.__Fprime = self.__X.transpose().dot(self.__F).dot(self.__X)

        # diagonalize Fock matrix
        try:
            self.__e, self.__Cprime = np.linalg.eigh(self.__Fprime)
        except np.linalg.LinAlgError:
            print('Error: eigenvalue convergence failed in iteration: %i' % niter)
            print('F:', self.__Fprime)
            print('H:', self.__H)
            print('J:', self.__J)
            print('XC:', self.__XC)
            print('P:', self.__P)
            raise np.linalg.LinAlgError
        
        # back-transform
        self.__C = self.__X.dot(self.__Cprime)
        
        # calculate total electronic energy
        M = self.__T + self.__V + 0.5 * self.__J
        P = self.__calculate_P()
        energy = np.einsum('ji,ij', P, M) + self.__Exc + self.__enuc
        
        # for the first few iterations, build a new density
        # matrix from the coefficients, else, resort to the DIIS
        # algorithm
        if niter <= SUBSPACE_START or not giis:
            P = self.__calculate_P()
            self.__P = (1.0 - mix) * self.__P + mix * P # use linear mixing

        # calculate DIIS coefficients
        e = (self.__F.dot(self.__P.dot(self.__S)) - \
             self.__S.dot(self.__P.dot(self.__F))).flatten()   # calculate error vector
        self.__enorm = np.linalg.norm(e)                       # store error vector norm
        self.__fmats_diis.append(self.__F)                     # add Fock matrix to list
        self.__pmat_diis.append(self.__P)                      # add density matrix to list
        self.__evs_diis.append(e)

        # prune size of the old Fock, density and error vector lists
        # only SUBSPACE_LENGTH iterations are used to guess the new
        # solution
        if len(self.__fmats_diis) > SUBSPACE_LENGTH:
            self.__fmats_diis = self.__fmats_diis[-SUBSPACE_LENGTH:]
            self.__pmat_diis = self.__pmat_diis[-SUBSPACE_LENGTH:]
            self.__evs_diis = self.__evs_diis[-SUBSPACE_LENGTH:]
        
        return energy
    
    def __setup(self):
        """
        Construct the bare classes and matrices necessary to start the
        self-consistent-field procedure
        """
        # construct basis functions and nuclei
        if issubclass(type(self.__basis), str): # if a basis set name is given
            self.__cgfs, self.__nuclei = self.__mol.build_basis(self.__basis)
        else: # either assume a list of CGFs objects is given
            self.__cgfs = self.__basis
            self.__nuclei = self.__mol.get_nuclei()
        
        # set number of electrons
        self.__nelec = np.sum([nucleus[1] for nucleus in self.__nuclei])

        # build molecular grid
        self.__molgrid = MolecularGrid([at for at in self.__mol],
                                       self.__cgfs, 
                                       nshells=self.__nshells, 
                                       nangpts=self.__nangpts,
                                       lmax=self.__lmax,
                                       fdpts=self.__fdpts,
                                       functional=self.__functional)
        self.__molgrid.initialize() # molecular grid uses late initialization

        # build one-electron matrices; because these matrices are Hermetian,
        # we only have to evaluate the upper half and then simply copy the
        # upper half to the lower half
        N = len(self.__cgfs)
        self.__S = np.zeros((N,N))
        self.__T = np.zeros_like(self.__S)
        self.__V = np.zeros_like(self.__S)
        for j in range(0,N):
            cgf1 = self.__cgfs[j]
            for i in range(j,N):
                cgf2 = self.__cgfs[i]
                self.__S[j,i] = self.__integrator.overlap(cgf1, cgf2)
                self.__T[j,i] = self.__integrator.kinetic(cgf1, cgf2)
                
                # in principle, the nuclear attraction could also be directly
                # obtained from the molecular grid, but analytical evaluation
                # is faster
                for k,nucleus in enumerate(self.__nuclei):
                    vjik = self.__integrator.nuclear(cgf1, cgf2, nucleus[0], nucleus[1])
                    self.__V[j,i] += vjik
                
                # copy upper triangle elements to lower triangle
                if i != j:
                    self.__S[i,j] = self.__S[j,i]
                    self.__T[i,j] = self.__T[j,i]
                    self.__V[i,j] = self.__V[j,i]
        
        # build single-electron matrix
        self.__H = self.__T + self.__V
        
        # diagonalize S and use it to construct the unitary transformation 
        # matrix that orthonormalizes the basis set
        s, U = np.linalg.eigh(self.__S)
        self.__X = U.dot(np.diag(1.0/np.sqrt(s)))
        
        # create empty matrices for the coulomb and exchange-correlation
        # energies
        self.__J = np.zeros_like(self.__S)
        self.__XC = np.zeros_like(self.__S)
        self.__Exc = 0.0
        
        # create zero density matrix
        self.__P = np.zeros_like(self.__S)
        
        # calculate nuclear repulsion
        self.__enuc = 0.0
        for i in range(0, len(self.__nuclei)):
            for j in range(i+1, len(self.__nuclei)):
                r = np.linalg.norm(self.__nuclei[i][0] - self.__nuclei[j][0])
                self.__enuc += self.__nuclei[i][1] * self.__nuclei[j][1] / r
        
        # build containers to store per-iteration data
        self.__energies = []
        self.__time_stats['iterations'] = []
        self.__fmats_diis = []
        self.__pmat_diis = []
        self.__evs_diis = []
        
    def __calculate_J(self):
        """
        Calculate the coulombic interaction matrix using the
        molecular grid
        """
        J, timestats = self.__molgrid.calculate_coulombic_matrix(True)
        self.calctimes['ulm_interpolation'].append(timestats['ulm_interpolation'])
        self.calctimes['build_hartree_field'].append(timestats['build_hartree_field'])
        self.calctimes['build_repulsion_matrix'].append(timestats['build_repulsion_matrix'])

        return J
    
    def __calculate_P(self):
        """
        Calculate density matrix from current coefficient matrix
        """
        N = len(self.__cgfs)
        P = np.zeros_like(self.__S)
        for i in range(0,N):
            for j in range(0,N):
                for k in range(0,int(self.__nelec//2)):
                    P[i,j] += 2.0 * self.__C[i,k] * self.__C[j,k]
                    
        return P
    
    def __calculate_XC(self):
        """
        Calculate the exchange-correlation matrix and the
        exchange-correlation energy
        """
        X, self.__Ex = self.__molgrid.calculate_exchange()
        C, self.__Ec = self.__molgrid.calculate_correlation()
        
        return X+C, self.__Ex + self.__Ec
               
    def __calculate_diis_coefficients(self, evs_diis):
        """
        Calculate the DIIS coefficients
        """
        B = np.zeros((len(evs_diis)+1, len(evs_diis)+1))
        B[-1,:] = -1
        B[:,-1] = -1
        B[-1,-1]=  0

        rhs = np.zeros((len(evs_diis)+1, 1))
        rhs[-1,-1] = -1

        for i in range(len(evs_diis)):
            for j in range(i+1):
                B[i,j] = np.dot(evs_diis[i].transpose(), evs_diis[j])
                B[j,i] = B[i,j]

        *diis_coeff, _ = np.linalg.solve(B,rhs)

        return diis_coeff

    def __extrapolate_fock_from_diis_coefficients(self, fmats_diis, diis_coeff):
        """
        Extrapolate the Fock matrix from the DIIS coefficients
        """
        norbs = fmats_diis[-1].shape[0]
        fguess = np.zeros((norbs,norbs))

        for i in range(len(fmats_diis)):
            fguess += fmats_diis[i]*diis_coeff[i]

        return fguess
        