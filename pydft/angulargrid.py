# -*- coding: utf-8 -*-

import os
import numpy as np
import pylebedev

class AngularGrid:
    def __init__(self):
        """
        Construct the class
        """
        
        # check if PyLebedev is available
        self.__leblib = pylebedev.PyLebedev()
        
        self.__lebedev_orders = [3, 5, 7, 9, 11, 13, 15, 17, 19, 
                                 21, 23, 25, 27, 29, 31, 35, 41, 47, 
                                 53, 59, 65, 71, 77, 83, 89, 95, 101, 
                                 107, 113, 119, 125, 131]

        # lebedev number of points
        self.__numpoints = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 
                            230, 266, 302, 350, 434, 590, 770, 974, 1202, 
                            1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 
                            4334, 4802, 5294, 5810]
        
        self.__coeff = dict()
    
    def get_coefficients(self, numpoints:int):
        """
        Get the coefficients for Lebedev quadrature on the unit sphere
        by number of points of that order
        """
        if numpoints in self.__numpoints:
            # load the data points from file if they do not yet exist
            if numpoints not in self.__coeff.keys():
                order = int(self.__lebedev_orders[self.__numpoints.index(numpoints)])
                self.__add_coefficients(order)
                
            # return datapoints
            return self.__coeff[str(numpoints)]
        else:
            raise Exception('There is no Lebedev order with %s number of points' % numpoints)
    
    def get_dataset_sizes(self):
        """
        Get the number of points on the unit sphere for each of the
        Lebedev orders
        """
        return [int(key) for key in self.__coeff.keys()]
                    
    def __add_coefficients(self, order):
        """
        Add coefficients from file and store these in a dictionary
        """
        points, weights = self.__leblib.get_points_and_weights(order)
        angles, weights = self.__leblib.get_points_and_weights(order, solid_angles = True)
        
        data = np.vstack([angles[:,0], angles[:,1], points[:,0], points[:,1], points[:,2], weights]).transpose()
        key = str(len(data))
        self.__coeff[key] = data