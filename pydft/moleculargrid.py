# -*- coding: utf-8 -*-

import numpy as np
from .atomicgrid import AtomicGrid
from . import bragg_slater
from .spherical_harmonics import spherical_harmonic
from pyqint import PyQInt
import pyqint as pq
import time
from .xcfunctionals import Functionals

class MolecularGrid:
    def __init__(self, 
                 atoms:list, 
                 cgfs:list[pq.cgf], 
                 nshells:int=32, 
                 nangpts:int=110, 
                 lmax:int=8,
                 functional:str='svwn5'):
        """
        Construct MolecularGrid

        Parameters
        ----------
        atoms : list
            list of atoms and their charges
        cgfs : list[pq.cgf]
            list of contracted Gaussian functions (basis set)
        nshells : int, optional
            number of radial shells, by default 32
        nangpts : int, optional
            number of angular sampling points per shell, by default 110
        lmax : int, optional
            maximum value for l for projection of spherical harmonics, by default 8
        functional : str, optional
            exchange-correlation functional, by default 'svwn5'
        """
        # keep track of build times
        self.construct_times = {}

        # set properties
        self.__atoms = atoms
        self.__nelec = np.sum([nuc[1] for nuc in self.__atoms])
        self.__lmax = lmax
        self.__nshells = nshells
        self.__nangpts = nangpts
        self.__basis = cgfs
        self.__functionals = Functionals(functional)
        self.__is_initialized = False
    
    def initialize(self):
        """
        Initialize the molecular grid
        
        This is a relatively expensive procedure and thus the user can
        delay the initialization of the molecular grid until.
        """
        if self.__is_initialized:
            return

        self.__build_molecular_grid()
        self.__build_amplitudes()
        self.__is_initialized = True

    def build_density(self, P:np.ndarray, normalize:bool=True):
        """
        Build the electron density from the density matrix

        Parameters
        ----------
        P : np.ndarray
            density matrix
        normalize : bool, optional
            whether to normalize the density matrix based on the total number
            of electrons, by default True
        """
        # this function requires the amplitudes and molecular grid
        # to be built
        self.initialize()


        # calculate densities at each grid point
        self.__densities = np.einsum('ijk,jl,ilk->ik', 
                                     self.__amplitudes, 
                                     P,
                                     self.__amplitudes)
        
        # also build the gradient of the density
        self.__gradients = 2.0 * np.einsum('ijk,jl,ilkm->ikm', 
                                           self.__amplitudes, 
                                           P,
                                           self.__ampgrads)
        
        # perform optional normalization
        if normalize:
            nelec = np.sum(self.__mgw * self.__densities.flatten())
            cn = self.__nelec / nelec  # normalization constant
            self.__densities *= cn
            self.__gradients *= cn
        
        # and place the densities back into the atomic grids
        for i,atgrid in enumerate(self.__atomgrids):
            atgrid.set_density(self.__densities[i,:])
            atgrid.set_gradient(self.__gradients[i,:,:])
            atgrid.build_hartree_potential()

    def get_becke_weights(self) -> np.ndarray:
        """
        Get the Becke coefficients

        Returns
        -------
        np.array
            :math:`(N \\times G \\times 3)` array with :math:`N` number of atoms 
            and :math:`G` number of grid points
        """
        return self.__mweights
    
    def get_grid_coordinates(self) -> list[np.array]:
        """
        Get the grid coordinates

        Returns
        -------
        list[np.array]
            list over atoms with for each atom a :math:`(G \\times 1)` array of
            electron densities where :math:`G` is the number of grid points 
            per atom
        """
        return self.__gridpoints
    
    def get_hartree_potential(self) -> np.array:
        """
        Get the Hartree potential at each grid cell

        Returns
        -------
        np.array
            :math:`(G \\times 1)` array with :math:`G` number of molecular
            grid points
        """
        return self.__ugpts
    
    def get_densities(self) -> list[np.array]:
        """
        Get the densities at the grid points
        
        Returns
        -------
        list[np.array]
            list over atoms with for each atom a :math:`(G \\times 1)` array of
            electron densities where :math:`G` is the number of grid points 
            per atom
        """
        return [atgrid.get_density() for atgrid in self.__atomgrids]
    
    def get_gradients(self) -> list[np.array]:
        """
        Get the gradient magnitude of the electron density field
        
        Returns
        -------
        list[np.array]
            list over atoms with for each atom a :math:`(G \\times 1)` array of
            electron gradient magnitudes where :math:`G` is the number of grid 
            points per atom
        """
        return [atgrid.get_gradient() for atgrid in self.__atomgrids]

    def get_rho_lm_atoms(self) -> np.array:
        """
        Get the electron density projected onto spherical harmonics per atom

        Returns
        -------
        np.array
            :math:`\\rho_{n,lm}` where :math:`n` loops over the atoms and the
            pairs :math:`lm` over the spherical harmonics
        """
        rho_lm = np.ndarray((len(self.__atoms), self.__nshells, (self.__lmax+1)**2))
        
        for i,at in enumerate(self.__atomgrids):
            rho_lm[i,:,:] = (at.get_rho_lm().T * at.get_radial_grid()**2).T
            
        return rho_lm

    def get_atomic_grid(self, i:int) -> AtomicGrid:
        """
        Get an atomic grid of atom :math:`i`

        Parameters
        ----------
        i : int
            atom index

        Returns
        -------
        AtomicGrid
            :class:`pydft.AtomicGrid` of atom :math:`i`
        """
        return self.__atomgrids[i]

    def count_electrons(self) -> float:
        """
        Count number of electrons in the system

        Returns
        -------
        float
            Total number of electrons
        """
        return np.sum([atgrid.count_electrons() for atgrid in self.__atomgrids])
    
    def count_electrons_from_rho_lm(self) -> float:
        """
        Count number of electrons from the spherical harmonic expansion; this
        function mainly acts to verify that the spherical harmonic expansion
        works correctly
        
        Returns
        -------
        float
            Total number of electrons

        Notes
        -----
        The Becke weights are already used when generating the density
        coefficients and hence for the back-transformation, we should *not*
        multiply with the Becke weights
        """
        # get spherical harmonic expansion coefficients per atomic grid
        ngpts = self.__nshells * self.__nangpts # gridpoints per atom
        rho_lm = np.ndarray((len(self.__atoms), self.__nshells, (self.__lmax+1)**2))
        
        for i,at in enumerate(self.__atomgrids):
            rho_lm[i,:,:] = at.get_rho_lm()
            
        # reshape the rho_lm array such that spherical harmonic coefficients
        # are set for each grid point and only store the coefficients for
        # each atomic grid associated with its own atom
        rho_lm_gpts = np.ndarray((len(self.__atoms), (self.__lmax+1)**2, ngpts))
        ylm = np.zeros_like(rho_lm_gpts)
        for i,at in enumerate(self.__atomgrids):
            for j in range(0, (self.__lmax+1)**2):
                rho_lm_gpts[i,j,:] = np.outer(rho_lm[i,:,j],np.ones(self.__nangpts)).flatten()
                ylm[i,j,:] = self.__ylmgpts[i,j,i*ngpts:(i+1)*ngpts]
    
        return np.einsum('ijk,ijk,ik', rho_lm_gpts,
                                       ylm,
                                       np.array([an.get_weights() for an in self.__atomgrids]))
    
    def calculate_dfa_nuclear_attraction_local(self) -> float:
        """
        Get density functional approximation for the nuclear attraction energy
        using only the local potential for each atomic grid and ignoring the
        influence of atomic attraction due to other atoms
        
        Returns
        -------
        float
            Approximate nuclear attraction energy (Hartree)

        Notes
        -----
        This function will of course not give an accurate result and is only
        meant for educational purposes
        """
        return np.sum([atgrid.get_dfa_nuclear_local() for atgrid in self.__atomgrids])
       
    def calculate_dfa_nuclear_attraction_full(self) -> np.array:
        """
        Get density functional approximation for the nuclear attraction energy

        Returns
        -------
        float
            Total nuclear attraction energy (Hartree)
        """
        return np.einsum('i,i,i', self.__npot, 
                                  np.array([an.get_density() for an in self.__atomgrids]).flatten(),
                                  self.__mgw)
    
    def calculate_dfa_coulomb(self) -> float:
        """
        Get density functional approximation for the coulombic repulsion 
        between electrons

        Returns
        -------
        float
            Total coulombic repulsion (Hartree)
        """
        # for each grid point, collect the spherical harmonic expansion coefficients
        # of the hartree potential for all the atoms
        ulmgpts = np.ndarray((len(self.__atoms), (self.__lmax+1)**2, len(self.__rgridpoints[0])))
        for i,at in enumerate(self.__atomgrids):
            ulmgpts[i,:,:] = at.calculate_interpolated_ulm(self.__rgridpoints[i])
        
        # for each grid point determine the Hartree potential from the 
        # spherical harmonic coefficients with respect to each atomic center
        self.__ugpts = np.einsum('ijk,ik,ijk->k', ulmgpts, self.__rigridpoints, self.__ylmgpts)
        
        return 0.5 * np.einsum('i,i,i', self.__ugpts, 
                                        np.array([an.get_density() for an in self.__atomgrids]).flatten(),
                                        self.__mgw)
    
    def calculate_dfa_coulomb_no_interpolation(self) -> float:
        """
        Calculate the Coulombic self-repulsion in the limit of the sum of atomic centers,
        without any interpolation of the contribution of the other atoms

        Returns
        -------
        float
            Approximate coulombic repulsion (Hartree)
        """
        return 0.5 * np.sum([at.calculate_coulomb_energy() for at in self.__atomgrids])

    def calculate_coulombic_matrix(self) -> np.array:
        r"""
        Build the coulomb matrix :math:`\mathbf{J}`

        Returns
        -------
        np.array
            Coulomb matrix :math:`\mathbf{J}`
        """
        # for each grid point, collect the spherical harmonic expansion coefficients
        # of the hartree potential for all the atoms
        ulmgpts = np.ndarray((len(self.__atoms), (self.__lmax+1)**2, len(self.__rgridpoints[0])))
        for i,at in enumerate(self.__atomgrids):
            ulmgpts[i,:,:] = at.calculate_interpolated_ulm(self.__rgridpoints[i])
        
        # for each grid point determine the Hartree potential from the 
        # spherical harmonic coefficients with respect to each atomic center
        self.__ugpts = np.einsum('ijk,ik,ijk->k', ulmgpts, self.__rigridpoints, self.__ylmgpts)
        
        # construct coulombic repulsion matrix by integrating the interaction
        # of the hartree potential with the basis function amplitudes
        N = len(self.__basis)
        J = np.zeros((N,N))
        for i in range(0, N):
            for j in range(i, N):
                J[i,j] = np.einsum('i,i,i,i', self.__ugpts, 
                                              self.__fullgrid_amplitudes[i,:],
                                              self.__fullgrid_amplitudes[j,:],
                                              self.__mgw)
                if i != j:
                    J[j,i] = J[i,j]
        
        return J
        
    def calculate_dfa_kinetic(self) -> float:
        """
        Get density functional approximation for the kinetic energy

        Returns
        -------
        float
            Total kinetic energy (Hartree)
        """
        return np.sum([atgrid.get_dfa_kinetic() for atgrid in self.__atomgrids])
    
    def calculate_dfa_exchange(self):
        """
        Get density functional approximation for the exchange energy

        Returns
        -------
        float
            Total exchange energy (Hartree)
        """
        return np.sum([atgrid.get_dfa_exchange() for atgrid in self.__atomgrids])

    def calculate_exchange(self) -> (np.ndarray,float):
        r"""
        Build the exchange matrix and give the exchange energy for the molecule

        Returns
        -------
        np.ndarray,float
            Exchange matrix :math:`\mathbf{X}` and exchange energy
        """
        
        dens = np.array([atgrid.get_density() for atgrid in self.__atomgrids]).flatten()
        if self.__functionals.is_gga():
            grad = np.array([atgrid.get_gradient_squared() for atgrid in self.__atomgrids]).flatten()
            fx, vfx = self.__functionals.calc_x(dens, grad)
        else:
            fx, vfx = self.__functionals.calc_x(dens)
        
        # calculate exchange energy
        ex = np.einsum('i,i,i', self.__mgw, fx, dens)
        
        # build matrix and pre-calculate weights
        X = np.zeros((len(self.__basis), len(self.__basis)))
        
        # exchange parameters
        for i in range(0, len(self.__basis)):
            for j in range(i, len(self.__basis)):
                X[i,j] = np.einsum('i,i,i,i', vfx, 
                                              self.__fullgrid_amplitudes[i,:],
                                              self.__fullgrid_amplitudes[j,:],
                                              self.__mgw)
                
                if i != j:
                    X[j,i] = X[i,j]
        
        return X, ex
    
    def calculate_correlation(self) -> (np.ndarray, float):
        """
        Build the correlation matrix and give the correlation energy for the molecule

        Returns
        -------
        np.ndarray,float
            Correlation matrix and exchange energy
        """
        dens = np.array([atgrid.get_density() for atgrid in self.__atomgrids]).flatten()
        if self.__functionals.is_gga():
            grad = np.array([atgrid.get_gradient_squared() for atgrid in self.__atomgrids]).flatten()
            fc, vfc = self.__functionals.calc_c(dens, grad)
        else:
            fc, vfc = self.__functionals.calc_c(dens)
        
        # calculate correlation energy
        ec = np.einsum('i,i,i', self.__mgw, fc, dens)
        
        # build matrix and pre-calculate weights
        C = np.zeros((len(self.__basis), len(self.__basis)))
        
        # exchange parameters
        for i in range(0, len(self.__basis)):
            for j in range(i, len(self.__basis)):
                C[i,j] = np.einsum('i,i,i,i', vfc, 
                                              self.__fullgrid_amplitudes[i,:],
                                              self.__fullgrid_amplitudes[j,:],
                                              self.__mgw)
                
                if i != j:
                    C[j,i] = C[i,j]
        
        return C, ec

    def get_density_at_points(self, 
                              spoints:np.ndarray, 
                              P:np.ndarray) -> np.array:
        """
        Calculate the electron density at the points as given by spoints
        using density matrix P

        Parameters
        ----------
        spoints : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points
        P : np.ndarray
            Density matrix  :math:`(K \\times K)` with :math:`K` the number of
            basis functions

        Returns
        -------
        np.array
            Electron density specified at grid points :math:`(N \\times 1)` 
            with :math:`N` the number of grid points
        """
        # build the amplitudes at the specified points
        amps = np.ndarray((len(self.__basis), len(spoints)))
        for i,cgf in enumerate(self.__basis):
            amps[i,:] = np.array([cgf.get_amp(p) for p in spoints])
        
        dens = np.einsum('ik,ij,jk->k', amps, P, amps)
        
        return dens
    
    def get_amplitude_at_points(self, 
                                spoints:np.ndarray, 
                                c:np.ndarray) -> np.ndarray:
        """
        Calculate the wave function amplitude at the points as given by 
        spoints using solution vector c

        Parameters
        ----------
        spoints : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points
        c : np.ndarray
            Coefficient vector :math:`\\vec{c}` of a single molecular orbital

        Returns
        -------
        np.array
            Electron density specified at grid points :math:`(N \\times 1)` 
            with :math:`N` the number of grid points
        """
        # build the amplitudes at the specified points
        amps = np.ndarray((len(self.__basis), len(spoints)))
        for i,cgf in enumerate(self.__basis):
            amps[i,:] = np.array([cgf.get_amp(p) for p in spoints])
        
        wfamp = np.einsum('ik,i->k', amps, c)
        
        return wfamp
    
    def get_gradient_at_points(self, 
                               spoints:np.ndarray, 
                               P:np.ndarray):
        """
        Calculate the gradient of the electron density at the points as
        given by spoints using density matrix :math:`\\mathbf{P}`

        Parameters
        ----------
        spoints : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points
        P : np.ndarray
            Density matrix  :math:`(K \\times K)` with :math:`K` the number of
            basis functions

        Returns
        -------
        np.array
            Electron density gradients at grid points :math:`(N \\times 3)` 
            with :math:`N` the number of grid points
        """
        # build the amplitudes at the specified points
        amps = np.ndarray((len(self.__basis), len(spoints)))
        grads = np.ndarray((len(self.__basis), len(spoints), 3))
        for i,cgf in enumerate(self.__basis):
            amps[i,:] = np.array([cgf.get_amp(p) for p in spoints])
            grads[i,:,:] = np.array([cgf.get_grad(p) for p in spoints])
        
        #
        # note that the 'ij' derivative component is equal
        # to the 'ji' derivative component
        #
        #t1 = np.einsum('ij,ikl,jk->kl', P, grads, amps)
        #t2 = np.einsum('ij,jkl,ik->kl', P, grads, amps)
        #return t1 + t2
        
        return 2.0 * np.einsum('ij,ikl,jk->kl', P, grads, amps)

    def calculate_coulomb_potential_at_points(self, 
                                              pts:np.ndarray) -> np.ndarray:
        r"""
        Collect the spherical harmonic expansion coefficients of the hartree 
        potential for all the atoms at the grid points

        Parameters
        ----------
        pts : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points

        Returns
        -------
        np.array
            Coulomb potential at grid points :math:`(\left[l_{\\textrm{max}}+1\\right]^{2} \\times N)` 
            with :math:`N` the number of grid points
        """
        # for each grid point, 
        ulmgpts = np.zeros(((self.__lmax+1)**2, len(pts)))
        for i,at in enumerate(self.__atomgrids):
            rgridpts = np.linalg.norm(pts - self.__atoms[i][0], axis=1)
            igpts = 1.0 / rgridpts
            val = at.calculate_interpolated_ulm(rgridpts)
            for j in range(0, len(val)):
                val[j,:] *= igpts
            ulmgpts += val
        
        return ulmgpts

    def get_exchange_potential_at_points(self, 
                                         pts:np.ndarray, 
                                         P:np.ndarray) -> np.ndarray:
        """
        Get the exchange potential at points pts

        Parameters
        ----------
        pts : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points
        P : np.ndarray
            Density matrix  :math:`(K \\times K)` with :math:`K` the number of
            basis functions

        Returns
        -------
        np.ndarray
            Exchange potential :math:`\\nu_{x}` at grid points :math:`(N \\times 1)` 
            with :math:`N` the number of grid points
        """
        dens = self.get_density_at_points(pts, P)
        if self.__functionals.is_gga():
            grad = np.linalg.norm(self.get_gradient_at_points(pts, P), axis=1)
            fx, vfx = self.__functionals.calc_x(dens, grad)
        else:
            fx, vfx = self.__functionals.calc_x(dens)
        
        return vfx
    
    def get_correlation_potential_at_points(self, 
                                            pts:np.ndarray, 
                                            P:np.ndarray) -> np.ndarray:
        """
        Get the correlation potential at points pts

        Parameters
        ----------
        pts : np.ndarray
            Array of grid points :math:`(N \\times 3)` with :math:`N` the number
            of grid points
        P : np.ndarray
            Density matrix  :math:`(K \\times K)` with :math:`K` the number of
            basis functions

        Returns
        -------
        np.ndarray
            Correlation potential :math:`\\nu_{c}` at grid points :math:`(N \\times 1)` 
            with :math:`N` the number of grid points
        """
        dens = self.get_density_at_points(pts, P)
        if self.__functionals.is_gga():
            grad = np.linalg.norm(self.get_gradient_at_points(pts, P), axis=1)
            fc, vfc = self.__functionals.calc_c(dens, grad)
        else:
            fc, vfc = self.__functionals.calc_c(dens)
        
        return vfc

    def get_spherical_harmonic_expansion_of_amplitude(self, 
                                                      c:np.ndarray, 
                                                      radial_factor=False) -> list[np.ndarray]:
        """
        Get the spherical harmonic expansion representation of a wavefunction
        amplitude using solution vector c

        Parameters
        ----------
        c : np.ndarray
            Wave function expansion coefficient
        radial_factor : bool, optional
            whether to multiply result by :math:`r^{2}`, by default False

        Returns
        -------
        list[np.ndarray]
            List of spherical harmonic expansions per atom stored as
            :math:`(N_{r} \\times N_{lm})) arrays where :math:`N_{r}` is the
            number of radial points and :math:`N_{lm}` the set of spherical
            harmonics used in the expansion.
        """
        she = []
        for at in self.__atomgrids:
            atgpts = at.get_full_grid()
            nr, na, _dummy = atgpts.shape
            amps = self.get_amplitude_at_points(atgpts.reshape(-1,3), c)
            she_at = at.perform_spherical_harmonic_expansion(amps.reshape((nr,na)))
            
            if radial_factor:
                r2 = at.get_radial_grid()**2
                she_at = (she_at.T * r2).T
            
            she.append(she_at)
        
        return np.array(she)

    def calculate_weights_at_points(self, 
                                    points:np.ndarray, 
                                    k:int=3) -> np.ndarray:
        """
        Custom function that (re-)calculates weights at given points

        Parameters
        ----------
        points : np.ndarray
            :math:`N \\times 3` array of grid points
        k : int, optional
            smoothing value, by default 3

        Returns
        -------
        np.ndarray
            :math:`N_{a} \\times N` array with :math:`N_{a}` the number of
            atoms and :math:`N` the number of grid points
        """
        mweights = np.zeros((len(self.__atoms), len(points)))
        smats = np.ndarray((len(self.__atoms), len(self.__atoms), len(points)))
        
        # loop over atoms
        for i in range(0, len(self.__atoms)):
            R1 = np.array(self.__atoms[i][0]) # position of atom 1
            rm1 = bragg_slater.BSRADII[self.__atoms[i][1]-1] # bragg-slater radius
            for j in range(i+1, len(self.__atoms)):
                R2 = np.array(self.__atoms[j][0]) # position of atom 2
                rm2 = bragg_slater.BSRADII[self.__atoms[j][1]-1] # bragg-slater radius
                Rij = np.linalg.norm(R1 - R2)  # distance between the two atoms
                chi = rm1 / rm2                # fraction of bragg-slater radii
                
                # calculate the eliptical coordinate mu at the point r on the
                # grid with respect to the two atoms
                muv = (np.linalg.norm(points - R1, axis=1) - \
                       np.linalg.norm(points - R2, axis=1)) / Rij
                
                # calculate the fuzzy cell coefficient between the two atoms i and j
                smats[i,j,:] = self.__sk(self.__vij(muv, chi), k)
                
                # automatically calculate the for j,i by taking one minus the result
                smats[j,i,:] = np.ones(len(smats[i,j,:])) - smats[i,j,:]
        
        # construct the cell functions for the points
        Pn = np.ones((len(self.__atoms), len(points)))
        for i in range(0, len(self.__atoms)):
            for j in range(0, len(self.__atoms)):
                if j != i:
                    Pn[i,:] = np.multiply(Pn[i,:], smats[i,j,:])
                    
        # calculate total of cell function for each atom at each grid point
        Pt = np.einsum('ij->j', Pn)
        
        # and normalize each point so that the sum equals unity
        for i in range(0, len(self.__atoms)):
            mweights[i,:] = np.divide(Pn[i,:], Pt) # cell function for each atom
            
        return mweights

    def __build_molecular_grid(self):
        """
        Build the molecular grid from the atomic grids
        """
        
        # build atomic grids
        st = time.time()
        self.__atomgrids = []
        for atom in self.__atoms:
            self.__atomgrids.append(AtomicGrid(atom, 
                                               self.__nshells, 
                                               self.__nangpts,
                                               self.__lmax)
                                    )
        self.construct_times['atomic_grids'] = time.time() - st
        
        # assign to each gridpoint in the atomic grid a weight in the molecular
        # grid
        st = time.time()
        self.__gridpoints = []
        self.__mweights = np.zeros((len(self.__atoms), 
                                    self.__nshells * int(self.__nangpts)))
        
        for g,atgrid in enumerate(self.__atomgrids):
            grid = np.array(atgrid.get_gridpoints())
            self.__gridpoints.append(grid)
            smats = np.ndarray((len(self.__atoms), len(self.__atoms), len(grid)))
            for i in range(0, len(self.__atoms)):
                R1 = np.array(self.__atoms[i][0]) # position of atom 1
                rm1 = bragg_slater.BSRADII[self.__atoms[i][1]-1] # bragg-slater radius
                for j in range(i+1, len(self.__atoms)):
                    R2 = np.array(self.__atoms[j][0]) # position of atom 2
                    rm2 = bragg_slater.BSRADII[self.__atoms[j][1]-1] # bragg-slater radius
                    Rij = np.linalg.norm(R1 - R2)  # distance between the two atoms
                    chi = rm1 / rm2                # fraction of bragg-slater radii
                    
                    # calculate the confocal elliptical coordinate mu at the point r on the
                    # grid with respect to the two atoms
                    muv = (np.linalg.norm(grid - R1, axis=1) - \
                           np.linalg.norm(grid - R2, axis=1)) / Rij
                    
                    # calculate the fuzzy cell coefficient between the two atoms i and j
                    smats[i,j,:] = self.__sk(self.__vij(muv, chi), 3)
                    
                    # automatically calculate the for j,i by taking one minus the result
                    smats[j,i,:] = np.ones(len(smats[i,j,:])) - smats[i,j,:]
            
            # for each atom, construct the cell function based on the product
            # of the fuzzy cells coefficient for that cell with each other cell
            Pn = np.ones((len(self.__atoms), len(grid)))
            for i in range(0, len(self.__atoms)):
                for j in range(0, len(self.__atoms)):
                    if j != i:
                        Pn[i,:] = np.multiply(Pn[i,:], smats[i,j,:])
            
            # calculate total of cell function for each atom at each grid point
            Pt = np.einsum('ij->j', Pn)
            
            # and normalize each point so that the sum equals unity
            for i in range(0, len(self.__atoms)):
                Pn[i,:] = np.divide(Pn[i,:], Pt) # cell function for each atom
                
            # store the result as a private variable
            self.__mweights[g,:] = Pn[g,:]
        
            # and store the molecular weight functions into the atomic grids
            atgrid.set_molecular_weights(self.__mweights[g])
        self.construct_times['fuzzy_cell_decomposition'] = time.time() - st
            
        # for each grid point in the molecular grid, store the distance to all the
        # nuclei in the grid
        st = time.time()
        self.__rgridpoints = np.zeros((len(self.__atoms), np.prod(self.__mweights.shape)))
        self.__xgridpoints = np.zeros_like(self.__rgridpoints)
        self.__ygridpoints = np.zeros_like(self.__rgridpoints)
        self.__zgridpoints = np.zeros_like(self.__rgridpoints)
        gpts = np.array(np.vstack(self.__gridpoints)) # global grid points
        for i,at in enumerate(self.__atoms):
            ap = at[0]
            # build relative coordinates
            self.__xgridpoints[i,:] = np.subtract(gpts[:,0], ap[0])
            self.__ygridpoints[i,:] = np.subtract(gpts[:,1], ap[1])
            self.__zgridpoints[i,:] = np.subtract(gpts[:,2], ap[2])
            
            # calculate radial distance between global gridpoint and atom i
            xy = np.add(np.power(self.__xgridpoints[i,:], 2),
                        np.power(self.__ygridpoints[i,:], 2))
            xyz = np.add(xy, np.power(self.__zgridpoints[i,:], 2))
            self.__rgridpoints[i,:] = np.sqrt(xyz)
        
        # build unit sphere angles
        self.__theta_gridpoints = np.arctan2(self.__ygridpoints,
                                             self.__xgridpoints)
        self.__phi_gridpoints = np.arccos(np.divide(self.__zgridpoints,
                                                    self.__rgridpoints))

        # calculate the nuclear potential at each grid point due to the other nuclei            
        self.__rigridpoints = np.divide(1.0, self.__rgridpoints)
        
        self.__npot = np.einsum('i,ij->j',
                                [-at[1] for at in self.__atoms], 
                                self.__rigridpoints)
        self.construct_times['nuclear_distance_and_potential'] = time.time() - st
        
        # calculate values for the spherical harmonics for each grid point
        # with respect to each atom as central point and for each value
        # of l,m (thus a rank-3 tensor)
        st = time.time()
        self.__ylmgpts = np.ndarray((len(self.__atoms), 
                                         (self.__lmax+1)**2, 
                                         np.prod(self.__mweights.shape)))
        for i,at in enumerate(self.__atoms):
            lmctr = 0
            for l in range(0, self.__lmax+1):
                for m in range(-l, l+1):
                    self.__ylmgpts[i,lmctr,:] = spherical_harmonic(l, m, \
                                                self.__theta_gridpoints[i,:], 
                                                self.__phi_gridpoints[i,:])
                    lmctr += 1
        self.construct_times['spherical_harmonics'] = time.time() - st
                    
        # collect complete weights for the full molecular grid
        weights = np.array([an.get_weights() for an in self.__atomgrids]).flatten()
        self.__mgw = np.multiply(weights, self.__mweights.flatten())

    def __build_amplitudes(self):
        """
        Precalculate the basis function amplitudes at the selected grid
        points
        """
        # calculate the amplitudes of the basis functions and store these
        # per atomic grid, per basis function and per points in the atomic
        # grid (i.e. in a rank-3 tensor)
        st = time.time()
        integrator = PyQInt()
        
        self.__amplitudes = np.ndarray((len(self.__atoms), len(self.__basis), len(self.__gridpoints[0])))
        for i,at in enumerate(self.__atoms):
            for j,cgf in enumerate(self.__basis):
                self.__amplitudes[i,j,:] = integrator.plot_wavefunction(np.array(self.__gridpoints[i]), [1], [cgf])
        
        # also build amplitude gradients
        self.__ampgrads = np.ndarray((len(self.__atoms), len(self.__basis), len(self.__gridpoints[0]), 3))
        for i,at in enumerate(self.__atoms):
            for j,cgf in enumerate(self.__basis):
                self.__ampgrads[i,j,:,:] = integrator.plot_gradient(np.array(self.__gridpoints[i]), [1], [cgf])

        # reassemble the amplitudes per atomic grid into one for the complete
        # molecular grid
        self.__fullgrid_amplitudes = np.ndarray((len(self.__basis), np.prod(self.__mweights.shape)))        
        for i,cgf in enumerate(self.__basis):
            self.__fullgrid_amplitudes[i,:] = np.hstack([self.__amplitudes[a,i,:] for a in range(len(self.__atoms))])
        self.construct_times['basis_function_amplitudes'] = time.time() - st

    def __step(self, mu):
        """
        Becke fuzzy grid cut-off function
        """
        if mu <= 0.0:
            return 1.0
        else:
            return 0.0

    def __vij(self, mu, chi):
        """
        Becke aij parameter for heteroatomic interactions
        """
        uij = (chi - 1.0) / (chi + 1.0)
        aij = uij / (uij**2 - 1.0)
        aij = min(aij, 0.5)
        aij = max(aij, -0.5)
        return mu + aij * (1.0 - mu**2)

    def __sk(self, mu, n):
        """
        Becke fuzzy grid cutoff profile
        """
        return 0.5 * (1.0 - self.__fuzzy(mu,n))

    def __fuzzy(self, mu, n):
        """
        Becke fuzzy grid iterative function for adaptive sharpening of
        the transition
        """
        if n == 0:
            return mu
        else:
            return self.__fuzzy(3.0/2.0 * mu - 0.5 * mu**3, n-1)