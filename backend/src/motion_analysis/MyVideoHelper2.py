#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 10:59:51 2017

@author: danifili
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import griddata

class MyVideoHelper2(object):
    
    def __init__(self, images):
        self.__images = images[:]
        
        ht = [-0.0378010678346322, 0.125047021427471, -0.267629124130554, 0.680287727944691]
        for i in range(4):
            ht.append(ht[3-i])
        self.__ht = ht
        
        gt = [-0.149035907898682, 0.204171130411112, -0.408622311811502, 1.69565453432942]
        for i in range(4):
            gt.append(-gt[3-i])
        self.__gt = gt
        
        hxy = [5.20468678432878e-5, 0.0369446223231693, 0.29028508996955, 0.536735013750202]
        for i in range(3):
            hxy.append(hxy[2-i])
        self.__hxy = hxy
        
        gxy = [0.00196396189835955, 0.112139472173319, 0.365478437146823, 1.59329214172416e-15]
        for i in range(3):
            gxy.append(-gxy[2-i])
        self.__gxy = gxy
        
        
    
    #creating properties
    def _width(self):
        return self.__images[0].width
    
    width = property(_width)
    
    def _height(self):
        return self.__images[0].height
    
    height = property(_height)
    
    def _duration(self):
        return len(self.__images)
    
    duration = property(_duration)
    
    def _images(self):
        return np.array(self.__images)
    
    images = property(_images)
    
    def __getitem__(self, coordinates):
        """
        Returns the brightness in a pixel given its coordinates
        
        params:
            -coordinates: A two dimensional tuple (x, y, t) of integers in which:
                        * 0 <= x < self.width and x = 0 represents left of image
                        * 0 <= y < self.height and y = 0 represents the top of image
                        * 0 <= t < self.duration
                
                        or a tuple of slices (x0:x1, y0: y1, t) in which:
                            * 0 <= x0 < x1 <= self.width
                            * 0 <= y0 < y1 <= self.height
                            * 0 <= t < self.duration

        returns:
            If coordinates is a tuple (x, y, t), an integer between 0 and 255 (inclusive) 
            which represents the brightness of pixel will be return.
            
            Otherwise, return a numpy array of shape (x1-x0, y1-y0) in which the (x, y)
            component corresponds a integer from 0 to 255 (inclusive) that represents the
            brightness at pixel (x0+x, y0+y) in the image at time t
        """
        
        # fail fast if preconditions are not met
        assert len(coordinates) == 3, "coordinates must be of length 3"
        
        x, y, t = coordinates
        
        assert type(x) is int or type(x) is slice, "x-coorndinate must an integer or a slice"
        assert type(y) is int or type(y) is slice, "y-coordinate must be and integer"
        assert type(x) is type(y), "x-coordinate and y-coordinate must be of the same type"
        
        
        assert type(t) is int, "time must be an integer"
        assert 0 <= t < self.duration, "time must be between 0 and self.duration"
        
        return self.__images[t][x,y]
    
        
    def _partial_ROI_helper(self, min_corner, max_corner, t, offset, partial_derivative):
        """
        Given the coordinates of the corners of a specific region of interest, a particular time and a shift of the image,
        and a function representing a partial derivative, it returns a matrix containing this partials on the specified
        region of interest.
        
        params:
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
            
            -max_corner: A tuple of size 2 representing the (x, y) coordinates
            of the lower-right corner of our region of interest.
            
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: represents the time. 
            
            -offset: A tuple of size 2 representing the shift of the first image when computing the partials.
            
            -partial_derivative: A function that takes two matrices of the same size as inputs, in which the first
            matrix represents the brightness at time t and the second matrix represents brightness at time t+1.
            
        returns:
            A matrix of floats where the (x, y) entry represents the partial derivative of B(x, y, t) specified 
            by the input partial_derivative function.
        
        """
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        sx, sy = offset
        
        assert 0 <= x_min < x_max < self.width
        assert 0 <= y_min < y_max < self.height
        assert 0 <= t < self.duration - 1
        
        
#        x_rel_A_min = x_min-sx-1
#        x_rel_A_max = x_max+2-sx
#        
#        y_rel_A_min = y_min-sy-1
#        y_rel_A_max = y_max+2-sy
#        
#        x_rel_B_min = x_min-1
#        x_rel_B_max = x_max+2
#        
#        y_rel_B_min = y_min-1
#        y_rel_B_max = y_max+2

        x_rel_min = x_min-2
        x_rel_max = x_max+3
        
        y_rel_min = y_min-2
        y_rel_max = y_max+3
        
        B = self.__images[(t+1)%8:]
        B.extend(self.__images[:(t+1)%8])
        B.reverse()
        
        B = [a[x_rel_min: x_rel_max, y_rel_min: y_rel_max] for a in B]
    
        return partial_derivative(B)

    
    def _Ex(self, min_corner, max_corner, t, offset = (0, 0)):
        """
        Given the coordinates of the corners of a specific region of interest, a particular time and a shift of the image,
        it returns the matrix containing the partial derivatives of brightness with respect to the x-coordinate
        
        params: 
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
            
            -max_corner: A tuple of size 2 representing the (x, y) coordinates
            of the lower-right corner of our region of interest.
            
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: represents the time
            
            -offset: A tuple of size 2 representing the shift of the first image when computing the partials
        
        returns:
            A matrix of floats where the (x, y) entry represents the partial derivative with respect to the
            x-coordinate of B(x,y,t) after shifting the first image
        """
        
        
#        kernel = np.array([[1,-1,0],
#                           [1,-1,0],
#                           [0,0,0]]).T
    
        kernel = np.array([[self.__gxy[x] * self.__hxy[y] for y in range(7)] for x in range(7)])
        
        dBdx = lambda A: sum(signal.convolve2d(A[i], kernel, mode='same') * self.__ht[7-i] for i in range(8))[2:-2, 2:-2]
    
        return self._partial_ROI_helper(min_corner, max_corner, t, offset, dBdx)
    
    def _Ey(self, min_corner, max_corner, t, offset = (0,0)):
        """
        Given the coordinates of the corners of a specific region of interest, a particular time and a shift of the image,
        it returns the matrix containing the partial derivatives of brightness with respect to the y-coordinate
        
        params: 
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
            
            -max_corner: A tuple of size 2 representing the (x, y) coordinates
            of the lower-right corner of our region of interest.
            
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: represents the time
            
            -offset: A tuple of size 2 representing the shift of the first image when computing the partials
        
        returns:
            A matrix of floats where the (x, y) entry represents the partial derivative with respect to the
            y-coordinate of B(x,y,t) after shifting the first image        
        """
#        kernel = np.array([[1,1,0],
#                           [-1,-1,0],
#                           [0,0,0]]).T
#    
#        kernel = kernel/2

        
        kernel = np.array([[self.__hxy[x] * self.__gxy[y] for y in range(7)] for x in range(7)])
    
        dBdy = lambda A: sum(signal.convolve2d(A[i], kernel, mode='same') * self.__ht[7-i] for i in range(8))[2:-2, 2:-2]
                             
        return self._partial_ROI_helper(min_corner, max_corner, t, offset, dBdy)
    
    def _Et(self, min_corner, max_corner, t, offset = (0,0)):
        """
        Given the coordinates of the corners of a specific region of interest, a particular time and a shift of the image,
        it returns the matrix containing the partial derivatives of brightness with respect to time
        
        params: 
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
            
            -max_corner: A tuple of size 2 representing the (x, y) coordinates
            of the lower-right corner of our region of interest.
            
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: represents the time
            
            -offset: A tuple of size 2 representing the shift of the first image when computing the partials
        
        returns:
            A matrix of floats where the (x, y) entry represents the partial derivative with respect to time t
            of B(x,y,t) after shifting the first image     
        """
#        kernel = np.array([[1/2,1/2,0],
#                           [1/2,1/2,0],
#                           [0,0,0]]).T
#        
#        kernel = kernel/2

        kernel = np.array([[self.__hxy[x] * self.__hxy[y] for y in range(7)] for x in range(7)])
        
        dBdt = lambda A: sum(signal.convolve2d(A[i], kernel, mode='same') * self.__gt[7-i] for i in range(8))[2:-2, 2:-2]
                            
        return self._partial_ROI_helper(min_corner, max_corner, t, offset, dBdt)
    
    def _get_average_displacements(self, displacements):
        """
        Calculate the local average displacement matrix given a matrix containing
        the displacements of a region of interest
        
        params:
            -displacements: A 3-dimensional numpy array of floats containing the displacements
            of a certain region of interest, in which the value at (x, y, 0) represents the displacements
            in the x-direction and (x, y, 1) in the y-direction
            
        returns:
            A matrix in which (x, y, d) is the average of (x+1,y,d), (x, y+1, d), (x-1, y, d), (x, y-1, d). If
            one of this points is out of bounds, its value will be replaced by the displacement of (x, y, d).
             
        """
        u, v = displacements[:, :, 0], displacements[:, :, 1]
        width, height = u.shape
        
        # Define kernel for convolution                                         
        kernel = np.array([[1/12,1/6,1/12],
                           [1/6 ,0  ,1/6],
                           [1/12,1/6,1/12]]) 

    
        # Perform 2D convolution with input data and kernel 
        u_avg = signal.convolve2d(u, kernel, boundary='symm', mode='same', fillvalue = 0)
        v_avg = signal.convolve2d(v, kernel, boundary='symm', mode='same', fillvalue = 0)
        
        averages = np.zeros((width, height, 2))
        averages[:, :, 0] = u_avg
        averages[:, :, 1] = v_avg
        
        return averages  
    
    def _interpolate(self, width, height, points, displacements):
        """
        Given the displacements of some pixels in a region of interest, it extends this displacements
        to the other pixels in the region of interest
        
        params:
            -width: positive integer; width of the region of interest
            
            -height: positive integer; height of the region of interest
            
            -points: a numpy array of tuples of two integers; it contains the pixels whose displacements are known
            
            -displacements: a numpy array of tuples of two floats; displacements[i] represents the displacements
                            in the x and y direction of points[i]
        
        return:
            a three dimensional numpy array containing the displacements in the region of interest after interpolation
        """

        displacements_x, displacements_y = displacements[:, 0], displacements[:, 1]
        
        grid_x, grid_y = np.mgrid[0:width, 0:height]
        
        final_displacements_x = griddata(points, displacements_x, (grid_x, grid_y), method="linear", fill_value=0)
        final_displacements_y = griddata(points, displacements_y, (grid_x, grid_y), method="linear", fill_value=0)
        
        final_displacements = np.zeros((width, height,2))
        final_displacements[:,:,0] = final_displacements_x
        final_displacements[:,:,1] = final_displacements_y
        
        return final_displacements
    

        
        
    
        
        
    
    
        
        