# -*- coding: utf-8 -*-

import numpy as np

class Functionals:
    """
    Class holding all exchange-correlation functionals
    
    Note that all the functions ending in _deriv are the derivative
    of the energy with respect to rho.
    
    To calculate the potential, nu, one has to calculate
    deltaf/deltarho = df/drho * rho + f
    """
    def __init__(self, functional = 'svwn5'):
        functionals = {
            'svwn5': [self.__slater, self.__vwn5, 
                      self.__slater_deriv, self.__vwn5_deriv, False],
            'pbe': [self.__pbe_x, self.__pbe_c, 
                    self.__pbe_x_deriv, self.__pbe_c_deriv, True]
        }
        
        if functional not in functionals.keys():
            raise Exception('Illegal XC-functional requested.')
        else:
            # store exchange and correlation functional
            self.__xf = functionals[functional][0]
            self.__cf = functionals[functional][1]
            self.__dxf = functionals[functional][2]
            self.__dcf = functionals[functional][3]
            self.__gga = functionals[functional][4]

    def is_gga(self):
        """
        Returns whether the loaded XC functional is of GGA type
        """
        return self.__gga

    def calc_x(self, rho, grad=None):
        """
        Calculate exchange using specified correlation potential
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            rho = self.__parse_dens(rho)
            if self.__gga:
                grad = self.__parse_dens(grad)
                ex = self.__xf(rho, grad)
                fx = self.__dxf(rho, grad)
            else:
                ex = self.__xf(rho)
                fx = self.__dxf(rho)
            
            return ex, (fx * rho + ex)
    
    def calc_c(self, rho, grad=None):
        """
        Calculate correlation using specified correlation potential
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            rho = self.__parse_dens(rho)
            if self.__gga:
                ex = self.__cf(rho, grad)
                fx = self.__dcf(rho, grad)
            else:
                ex = self.__cf(rho)
                fx = self.__dcf(rho)
            
            ex = np.nan_to_num(ex)
            fx = np.nan_to_num(fx)
            
            return ex, (fx * rho + ex)

    def __parse_dens(self, dens):
        return np.maximum(dens, 1e-12)

    def __slater(self, dens):
        """
        Slater exchange functional
        """
        return -(3/4) * (3 / np.pi)**(1/3) * np.power(dens, 1./3.)
        
    def __slater_deriv(self, dens):
        """
        Slater exchange functional derivative towards rho
        """
        return -(3 / np.pi)**(1/3) / 4 / dens**(2/3)
        
    def __vwn5(self, dens):
        # Note that equation E.27 is incorrect in the book of Parr and Yang
        # but is correct in the original paper of Vosko, Wilk and Nusair
        # https://cdnsciencepub.com/doi/pdf/10.1139/p80-159
        A = 0.0621814
        x0 = -0.409286
        b = 13.0720
        c = 42.7198
        
        rs = (3. / 4. / np.pi / dens)**(1./3.)
        
        x = rs**(1/2)
        X = x**2 + b * x + c
        X0 = x0**2 + b * x0 + c
        Q = (4 * c - b**2)**(1/2)
        atan = np.arctan(Q / (2*x+b))
        
        return A/2 * (np.log(x**2/X) + 2*b/Q * atan - b*x0/X0 * (np.log((x-x0)**2 / X) + 2 * (b + 2*x0) / Q * atan))

    def __vwn5_deriv(self, dens):
        """
        Derivative of the VWN5 correlation functional towards the density
        """
        A = 0.0621814
        x0 = -0.409286
        b = 13.0720
        c = 42.7198
        X0 = x0**2 + b * x0 + c
        Q = (4 * c - b**2)**(1/2)

        rs = (3. / 4. / np.pi / dens)**(1./3.)
        x = rs**(1/2)

        T1 = (2*c + b*x)/(x*(c + x*(b + x)))
        T2 = (4*b)/(Q**2 + (b + 2*x)**2)
        T3 = (b*x0*(-((b + 2*x)/(c + x*(b + x))) + 2/(x - x0) - (4*(b + 2*x0))/(Q**2 + (b + 2*x)**2))) / X0
        
        dfdx = A/2 * (T1 - T2 - T3)
        dxdrho = -((1/dens)**(7/6)/(2*2**(1/3)*3**(5/6)*np.pi**(1/6)))
        return dfdx * dxdrho
    
    def __pbe_x(self, rho, gamma):
        """
        PBE exchange functional for spin-unpolarized density
        """
        R = 0.804
        mu = 0.2195149727645171
        
        c1 = 0.738558766382
        c2 = 0.0192920212964
        c3 = 0.0261211729852
        
        return rho**(1/3) * (-c1 - c2 * gamma * mu * R \
                             / (c3 * gamma * mu + R * rho**(8/3)))

    def __pbe_x_deriv(self, rho, gamma):
        """
        PBE exchange derivative for spin-unpolarized density using reparametrization
        """
        R = 0.804
        mu = 0.2195149727645171

        c1 = 0.246186255461
        c2 = 18.8495559215
        c3 = 65.9734457254
        c4 = 360.809905669
        c5 = 38.2831200025
        
        t1 = gamma**2 * mu**2 * (-c1 - c1 * R)
        t2 = gamma * mu * R * (-c2 + c3 * R) * rho**(8/3)
        t3 = -c4 * R**2 * rho**(16/3)
        t4 = (rho**(2/3) * (gamma * mu + c5 * R * rho**(8/3))**2)
        
        return (t1 + t2 + t3) / t4

    def __pbe_x_deriv_numerical(self, rho, gamma):
        """
        PBE Exchange derivative using numerical approximation
        """
        # use finite difference discretization to approximate result
        dx = 1e-5
        return (self.pbe_x_deriv_simplified(rho + dx, gamma) - self.pbe_x_deriv_simplified(rho - dx, gamma)) / (2. * dx)

    def __pbe_c(self, rho, gamma):
        """
        Reparametrization of the PBE correlation functional for spin-unpolarized density
        """
        c1 = -0.0621814 
        c2 = 0.008243319792565314
        c3 = 0.3720033616668558
        c4 = 0.13838902240433937
        c5 = 0.049771773124349425
        c6 = 0.011795838478103926
        c7 = 0.031090690869654897
        c8 = 7.341562356353513
        c9 = 0.13621078885675922
        c10 = 0.3720033616668558
        c11 = 0.049771773124349425
        c12 = 2.0000005873362636
        c13 = 0.2651378776729259
        
        irho = 1/rho
        irho12 = np.sqrt(irho)
        irho16 = (irho)**(1/6)
        irho13 = irho16**2
        irho23 = irho13**2
        
        c3irho16 = c3*irho16
        c6irho23 = c6*irho23
        c13irho13 = c13*irho13
        c4irho13 = c4*irho13
        c5irho12 = c5*irho12
        c9gamma = c9 * gamma
        
        t1 = (c1 - c2*irho13)
        t2 = 1 + 1/(c3irho16 + c4irho13 + c5irho12 + c6irho23)
        t3 = (-1. + 1.*(1 + 1/(c3irho16 + c4irho13 + c5irho12 + c6irho23))**(c12 + c13irho13))
        t4 = (-1. + 1.*(1 + 1/(c10*irho16 + c4irho13 + c11*irho12 + c6irho23))**(c12 + c13irho13))
        t5 = 1. + 1 / ((c8*rho**(7/3))/gamma + (c9gamma)/(t4*(c9gamma + t3*rho**(7/3))))

        return t1*np.log(t2) + c7*np.log(t5)

    def __pbe_c_deriv(self, rho, gamma):
        """
        PBE correlation derivative using numerical approximation
        """
        # use finite difference discretization to approximate result
        dx = 1e-5
        return (self.__pbe_c(rho + dx, gamma) - self.__pbe_c(rho - dx, gamma)) / (2. * dx)