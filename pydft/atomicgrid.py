# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import CubicSpline
from . import bragg_slater
from .angulargrid import AngularGrid
from .spherical_harmonics import spherical_harmonic

class AtomicGrid:
    def __init__(self, at, nshells:int=32, nangpts:int=110, lmax:int=8):
        """
        Initialize the class
        
        at:      atom
        nshells: number of radial shells
        nangpts: number of angular points
        """
        self.__atom = at
        self.__cp = at[0]   # center of the atomic grid (= position of the atom)
        self.__aidx = at[1] # element id
        self.__nshells = nshells
        
        # set coefficients and parameters
        self.__set_bragg_slater_radius()
        self.__build_chebychev_grid()
        self.__lebedev_coeff = self.__load_lebedev_coefficients(nangpts)
        
        # store the positions of the angular points, their weights for the
        # Lebedev integration and the value of the spherical harmonics at
        # these points
        self.__usangles = [np.array([p[0], p[1]]) for p in self.__lebedev_coeff]
        self.__angpts = [np.array([p[2], p[3], p[4]]) for p in self.__lebedev_coeff]
        self.__wangpts = np.array([p[5] for p in self.__lebedev_coeff])
        self.__lmax = lmax
        self.__build_ylm()
        self.__build_weight_grid()
        self.__mweights = np.ones_like(self.__wgrid)
    
    def get_charge(self):
        """
        Get the atomic charge
        """
        return self.__aidx
    
    def get_radial_grid(self):
        """
        Return radial grid
        """
        return self.__rr
    
    def get_full_grid(self):
        """
        Get all the grid points as a list of positions (Nx3) matrix
        """
        return np.multiply.outer(self.__rr, self.__angpts)
    
    def set_density(self, density):
        """
        Set the density at the grid points
        """
        self.__edens = density.reshape((len(self.__rr), len(self.__angpts)))
        
    def set_gradient(self, gradient):
        """
        Set the density at the grid points
        """
        self.__gradient = gradient.reshape((len(self.__rr), len(self.__angpts), 3))
    
    def get_local_hartree_potential(self):
        """
        Returns the local Hartree potential for only this atomic grid and
        not considering any interactions of the other atomic grids
        
        Note:
        * This function should only be used after htpot has been generated
          using either calculate_coulomb_energy() or variants of this function
        * This function should only be used for educational purposes
        """
        return self.__htpot.flatten()
    
    def build_hartree_potential(self):
        """
        Construct the Hartree potential at the grid points including cubic
        interpolation to obtain the Hartree potential expansion coefficients
        at other positions
        """
        # for each radial grid point, construct the set of spherical harmonic
        # coefficients that describes the density in that radial shell
        self.__calculate_density_coefficients()
        
        # from the radial coefficients (rho_lm), construct the spherical
        # harmonic coefficients that describe the Hartree potential (ulm)
        # in that radial shell
        self.__calculate_hartree_potential_coefficients()
        
        # build cubic interpolation
        self.__ulm_splines = []
        rr = np.flip(self.__rr)
        for i in range((self.__lmax+1)**2):
            self.__ulm_splines.append(CubicSpline(rr, np.flip(self.__ulm[:,i])))

    def get_ulm(self):
        """
        Get Hartree potential coefficients
        """
        return self.__ulm

    def get_rho_lm(self):
        """
        Get electron density coefficients
        """
        return self.__rho_lm
    
    def get_ylm(self):
        """
        Get values for the spherical harmonics
        """
        return self.__ylm

    def calculate_interpolated_ulm(self, rr):
        """
        Calculate interpolated Hartree potential expansion coefficients
        for all points r in the array rr
        """
        ulmintp = np.zeros((len(self.__ulm_splines), len(rr)))
                
        # loop over lm pairs and interpolate value for Ulm at points rr
        for i in range((self.__lmax+1)**2):
            ulmintp[i,:] = self.__ulm_splines[i](rr)

        return ulmintp

    def set_molecular_weights(self, mweights):
        """
        Set the molecular weights (Becke weights) at the grid points
        """
        self.__mweights = mweights.reshape((len(self.__rr), len(self.__angpts)))

    def get_gridpoints(self):
        """
        Return list of grid points
        """
        return self.__gridpoints
    
    def get_weights(self):
        """
        Get the weights on the atomic grid for the quadrature
        """
        return self.__wgrid.flatten()
    
    def get_weights_angular_points(self):
        """
        Get the weights for the angular points only
        """
        return self.__wangpts
    
    def get_becke_weights(self):
        """
        Return the fuzzy cell weights (Becke weights)
        """
        return self.__mweights

    def get_density(self):
        """
        Get the number of electrons for this atomic grid
        """
        return self.__edens.flatten()
    
    def get_gradient(self):
        """
        Get the gradient of the electron density for this atomic grid
        """
        return self.__gradient.reshape((len(self.__rr) * len(self.__angpts), 3))
    
    def get_gradient_squared(self):
        """
        Get the modulus squared of the gradient for this atomic grid
        """
        grad = self.get_gradient()
        #return np.linalg.norm(grad, axis=1)
        return np.einsum('ij,ij->i', grad, grad)

    def count_electrons(self):
        """
        Sum over the electron density to obtain the number of electrons for
        this atomic cell
        """
        return np.einsum('ij,ij,ij', self.__edens, self.__wgrid, self.__mweights)
    
    def get_dfa_exchange(self):
        """
        Get density functional approximation of the exchange energy for this atomic cell
        """
        alpha = 2.0 / 3.0
        fac = -2.25 * alpha * np.power(0.75 / np.pi, 1.0 / 3.0)
        return fac * np.einsum('ij,ij,ij', np.power(self.__edens, 4/3), self.__wgrid, self.__mweights)
    
    def get_dfa_kinetic(self):
        """
        Get density functional approximation of the kinetic energy for this atomic cell
        """
        Ckin = 3.0 / 40.0 * (3 / np.pi)**(2/3) * (2 * np.pi)**2
        return Ckin * np.einsum('ij,ij,ij', np.power(self.__edens, 5/3), self.__wgrid, self.__mweights)

    def get_dfa_nuclear_local(self):
        """
        Get density functional approximation of the kinetic energy for this atomic cell
        """
        npot = -self.__aidx / self.__rr # build nuclear attraction potential
        return np.einsum('i,ij,ij,ij', npot, self.__edens, self.__wgrid, self.__mweights)

    def perform_spherical_harmonic_expansion(self, f):
        """
        Calculate the expansion coefficients using a basis set composed of
        spherical harmonics given a atomic orbital
        
        The expansion is performed over the function f which needs to be a matrix
        of Nk x Nj elements where Nk is the number of radial points and Nj is
        the number of angular points.
        """
        return np.einsum('ij,kj,j,kj->ki', 
                         self.__ylm, 
                         f, 
                         self.__wangpts,
                         self.__mweights) * 4.0 * np.pi

    def get_bragg_slater_radius(self) -> float:
        """
        Get the Bragg-Slater radius
        
        Returns
        -------
        float
            Bragg-Slater radius
        """
        return self.__rm

    def __calculate_density_coefficients(self):
        """
        Calculate the expansion coefficients using a basis set composed of
        spherical harmonics given a atomic orbital
        
        Note that we explicitly multiply with the Becke weights so that we
        only integrate over the electron density assigned for this fuzzy cell.
        That means that in any back-transformation, we should explicitly not
        multiply with the Becke weights
        """
        self.__rho_lm = self.perform_spherical_harmonic_expansion(self.__edens)
    
    def __calculate_hartree_potential_coefficients(self):
        """
        Calculate the Hartree potential using the finite difference approximation
        """
        # construct base matrix
        Mbase = self.__build_finite_difference_matrix(self.__rr, 0.5 * self.__rm)
        self.__ulm = np.zeros_like(self.__rho_lm)

        lmctr = 0 # indexing counter
        for l in range(0, self.__lmax+1):
            for m in range(-l, l+1):
                M = Mbase.copy() # create copy of base matrix
                N = M.shape[0]   # grid size + 2
                b = np.zeros(N)
                
                for i in range(0, N):
                    # set first element for first lm value
                    if i == 0:
                        if l == 0:
                            b[i] = np.sqrt(4. * np.pi) * self.count_electrons()
                        continue
                    
                    # skip last element
                    if i == N-1:
                        continue
                
                    M[i,i] -= l * (l+1) / self.__rr[i-1]**2
                    b[i] = -4.0 * np.pi * self.__rr[i-1] * self.__rho_lm[i-1, lmctr]
                
                self.__ulm[:,lmctr] = np.linalg.solve(M,b)[1:-1]

                # increment lm counter                    
                lmctr += 1
    
    def calculate_nuclear_attraction(self):
        """
        Calculate the nuclear attraction given the electron density
        """
        npot = -self.__aidx / self.__rr # build nuclear attraction potential
        return np.einsum('ij,i,ij', self.__wgrid, npot, self.__edens)
    
    def calculate_coulomb_energy(self):
        """
        Calculate the electron-electron repulsion for this atomic grid
        
        This function is only for educational purposes, for any evaluation
        of the coulomb energy, extensive interpolation over all the atomic
        grids is necessary
        """
        # build Hartree potential
        pot = 1.0 / self.__rr
        self.__htpot = np.einsum('ij,jk,i,ik->ik', self.__ulm, self.__ylm, pot, self.__edens)
                    
        return np.einsum('ij,ij', self.__htpot, self.__wgrid)
    
    def calculate_coulomb_energy_interpolation(self):
        """
        Calculate the electron-electron repulsion for this atomic grid
        using cubic spline interpolation
        
        This function is only for educational purposes, for any evaluation
        of the coulomb energy, extensive interpolation over all the atomic
        grids is necessary
        """
        # build Hartree potential
        pot = 1.0 / self.__rr
        ulmintp = self.calculate_interpolated_ulm(self.__rr).transpose()
        self.__htpot = np.einsum('ij,jk,i,ik->ik', ulmintp, self.__ylm, pot, self.__edens)
                    
        return np.einsum('ij,ij', self.__htpot, self.__wgrid)
    
    def __build_weight_grid(self):
        """
        Calculates the weight factors on the grid and also builds the list
        of coordinates of the grid points
        """
        radial_weights = np.array([wgc * r**2 for wgc,r in zip(self.__wgcs, self.__rr)])
        self.__wgrid = np.outer(4.0 * np.pi * radial_weights, self.__wangpts)
        
        self.__gridpoints = []
        for rc,r in enumerate(self.__rr):
            [self.__gridpoints.append(r * angpt + self.__cp) for angpt in self.__angpts]
    
    def __load_lebedev_coefficients(self, order:int):
        """
        Load Lebedev coefficients from data file and store these as a class
        variable in a dictionary. Each set of angular points is stored as
        a Nx4 matrix wherein the first three per row correspond to the
        position on the unit sphere and the last row to the weight as used
        in the Lebedev integration.
        """
        ag = AngularGrid()
        return ag.get_coefficients(order)
    
    def __build_chebychev_grid(self):
        """
        Build Chebychev grid of concentric spheres
        """
        f = np.pi / float(self.__nshells+1)
        z = np.array(range(1, self.__nshells+1))
        x = np.cos(f * z)
        rm = self.__rm * 0.5
        
        self.__rr = rm * (1.0 + x) / (1.0 - x)
        self.__wgcs = f * np.sin(f * z)**2 / np.sqrt(1.0 - x**2) * 2.0 * rm / (1.0 - x)**2
    
    def __build_ylm(self):
        """
        Calculate the value for spherical harmonics at the angular points
        """
        self.__ylm = np.zeros(((self.__lmax+1)**2, len(self.__angpts)))
        lmctr = 0
        for l in range(0, self.__lmax+1):
            for m in range(-l, l+1):
                ylm = [spherical_harmonic(l, m, p[0], p[1]) for p in self.__usangles]
                self.__ylm[lmctr,:] = np.array(ylm)
                lmctr += 1
    
    def __build_finite_difference_matrix(self, r, rm):
        """
        Construct the finite difference matrix to solve for Ulm
        
        r : vector of distances
        rm: half of the Bragg-Slater radius
        """
        N = len(r)
        A = np.zeros((N+2, N+2))
        h = 1.0 / float(N+1)
        
        for i in range(0, N+2):
            c1 = 0.0
            c2 = 0.0
            
            if i > 0 and i < N+1:
                c1 = self.__dzdrsq(r[i-1], rm)
                c2 = self.__d2zdr2(r[i-1], rm)
            
            if i == 0:
                A[0,0] = 1.0
                continue
            
            if i == 1:
                c1 /= 12.0 * h * h
                c2 /= 12.0 * h
                A[i,0] = 11.0 * c1 -3.0 * c2
                A[i,1] = -20.0 * c1 - 10.0 * c2
                A[i,2] = 6.0 * c1 + 18.0 * c2
                A[i,3] = 4.0 * c1 - 6.0 * c2
                A[i,4] = -1.0 * c1 + 1.0 * c2
                continue
            
            if i == 2:
                c1 /= 12.0 * h * h
                c2 /= 60.0 * h
                A[i,0] = -1.0 * c1 + 3.0 * c2
                A[i,1] = 16.0 * c1 - 30.0 * c2
                A[i,2] = -30.0 * c1 - 20.0 * c2
                A[i,3] = 16.0 * c1 + 60.0 * c2
                A[i,4] = -1.0 * c1 - 15.0 * c2
                A[i,5] = 0.0 * c1 + 2.0 * c2
                continue
            
            if i == N-1:
                c1 /= 12.0 * h * h
                c2 /= 60.0 * h
                A[i,N-4] = 0.0 * c1 - 2.0 * c2
                A[i,N-3] = -1.0 * c1 + 15.0 * c2
                A[i,N-2] = 16.0 * c1 - 60.0 * c2
                A[i,N-1] = -30.0 * c1 + 20.0 * c2
                A[i,N] = 16.0 * c1 + 30.0 * c2
                A[i,N+1] = -1.0 * c1 - 3.0 * c2
                continue
            
            if i == N:
                c1 /= 12.0 * h * h
                c2 /= 12.0 * h
                A[i,N-3] = -1.0 * c1 - 1.0 * c2
                A[i,N-2] = 4.0 * c1 + 6.0 * c2
                A[i,N-1] = 6.0 * c1 - 18.0 * c2
                A[i,N] = -20.0 * c1 + 10.0 * c2
                A[i,N+1] = 11.0 * c1 + 3.0 * c2
                continue
            
            if i == N+1:
                A[i,i] = 1.0
                continue
            
            c1 /= 180.0 * h * h 
            c2 /= 60.0 * h
            A[i,i-3] = 2.0 * c1 - 1.0 * c2
            A[i,i-2] = -27.0 * c1 + 9.0 * c2
            A[i,i-1] = 270.0 * c1 - 45.0 * c2
            A[i,i]   = -490.0 * c1
            A[i,i+1] = 270.0 * c1 + 45.0 * c2
            A[i,i+2] = -27.0 * c1 - 9.0 * c2
            A[i,i+3] = 2.0 * c1 + 1.0 * c2
            
        return A

    def __set_bragg_slater_radius(self):
        """
        Set the Bragg-Slater radius for the atom
        """
        # note that the bragg_slater.BSRADII list starts with hydrogen at
        # index 0 whereas aidx variable sets Hydrogen at index 1, Helium at
        # index 2, and so on
        self.__rm = bragg_slater.BSRADII[self.__aidx-1]

    def __d2zdr2(self, r, rm):
        """
        Calculate second derivative of z-grid towards (regular) r-grid
        """
        nom = rm * rm * (rm + 3.0 * r)
        denom = 2.0 * np.pi * (rm * r)**(3/2) * (rm + r) ** 2
        return nom / denom

    def __dzdrsq(self, r, rm):
        """
        Calculate square of first derivative of z-grid towards (regular) r-grid
        """
        return rm / (np.pi * np.pi * r * (rm + r) * (rm + r))