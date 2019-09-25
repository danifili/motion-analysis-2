#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:16:20 2017

@author: danifili
"""
import numpy as np
from scipy import signal
from src.motion_analysis.DenseOpticalFlow import AbstractDenseOpticalFlow
from timeit import default_timer as timer
from tqdm import tqdm

class HornShunck(AbstractDenseOpticalFlow):
    """
    Represents a video of a sequence of images. It methods for estimating the optical flow
    are implemented with the Horn-Shunck method
    """
    def __init__(self, images, Ex=None, Ey=None, Et=None):
        """
        Creates a MyVideo3 object from a sequence of two or more images of
        the same size
        
        params:
          -images: a list of two or more MyImage objects of the same size
        """
        
        AbstractDenseOpticalFlow.__init__(self, images) 
        
        #precompute matrices of gradients
        min_corner = (0, 0)
        max_corner = (self.width-1, self.height-1)
        
        image_shape = (self.width, self.height, self.duration-1)
        
        if Ex is None:
            self.__Ex = np.zeros(image_shape)
        if Ey is None:
            self.__Ey = np.zeros(image_shape)
        if Et is None:
            self.__Et = np.zeros(image_shape)
        
        for t in range(self.duration-1):
            if Ex is None:
                self.__Ex[:, :, t] = self._Ex(min_corner, max_corner, t)
            if Ey is None:
                self.__Ey[:, :, t] = self._Ey(min_corner, max_corner, t)
            if Et is None:
                self.__Et[:, :, t] = self._Et(min_corner, max_corner, t)
        
        if Ex is not None:
            self.__Ex = np.array(Ex)
        if Ey is not None:
            self.__Ey = np.array(Ey)
        if Et is not None:
            self.__Et = np.array(Et)
        
        
    
    def get_optical_flow_ROI(self, min_corner, max_corner, t, initial_displacements = None, max_iterations = 10000,
                             smoothness = 0.1, input_mask = None):
        
        #get corners coordinates
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        #check preconditions are met
        assert 0 <= x_min < x_max < self.width
        assert 0 <= y_min < y_max < self.height
        assert 0 <= t < self.duration -1
        
        #get region of interest width and height
        roi_width = x_max+1 - x_min
        roi_height = y_max+1 - y_min
        
        #initialize numpy arrays containing  current and previous displacements
        current_displacements = np.zeros((roi_width, roi_height, 2))
        
        if input_mask is None:
            mask = np.zeros((roi_width, roi_height, 2))
        else:
            mask = np.array(input_mask)
            
        anti_mask = np.ones((roi_width, roi_height, 2)) - mask
        
        if initial_displacements is None:
            previous_displacements = np.zeros((roi_width, roi_height, 2))
        else:
            previous_displacements = np.array(initial_displacements)
        
        if max_iterations <= 0:
            return previous_displacements
            
        weighted_initial_displacements = np.multiply(mask, previous_displacements)
        
        #compute the partial derivatives with respect to space and time
        Ex = self.__Ex[x_min:x_max+1, y_min:y_max+1, t]
        Ey = self.__Ey[x_min:x_max+1, y_min:y_max+1, t]
        Et = self.__Et[x_min:x_max+1, y_min:y_max+1, t]
        
        den = smoothness ** 2 * np.ones((roi_width, roi_height)) + np.multiply(Ex, Ex) + np.multiply(Ey, Ey)
        a = np.divide(Ex, den)
        b = np.divide(Ey, den)
        c = np.divide(Et, den)
        
        start = timer()
        
        for each in tqdm(range(max_iterations)):
            #get update the displacements based on the previous displacements
            current_displacements = self._update_displacements(previous_displacements, a, b, c)
            current_displacements = np.multiply(anti_mask, current_displacements) + weighted_initial_displacements
            
            #swap contents
            current_displacements, previous_displacements = previous_displacements, current_displacements
        
        end = timer()
        print (end-start)
        #return matrix with displacements
        return previous_displacements


    def _update_displacements(self, displacements, a, b, c):
        """
        Update the pixel displacements given a matrix of displacements
        and matrices containing the gradients
        
        params:
            displacements: A 3-dimensional numpy array of floats containing the current displacements
            of a certain region of interest, in which the value at (x, y, 0) represents the current displacement 
            in the x-direction and (x, y, 1) in the y-direction
            
            Ex: the matrix of the partial derivatives of brightness with respect to the x-coordinate
            Ey: the matrix of the partial derivatives of brightness with respect to the y-coordinate
            Et: the matrix of the partial derivatives of brightness with respect to time
            
            -smoothness: positive flow used to control the smoothness of the optical flow. The closer to 0 
            the smoother

        returns:
            A matrix containing the new displacements
        """
        #get average displacements from neighbors
        averages = self._get_average_displacements(displacements)
        
        u_average, v_average = averages[:, :, 0], averages[:, :, 1]
        width, height = u_average.shape
                       
        error = np.multiply(a, u_average) + np.multiply(b, v_average) + c
        
        return averages - error[:,:,np.newaxis]


        
                 
    
