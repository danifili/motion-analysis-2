#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 11:25:52 2017

@author: danifili
"""
import numpy as np
from scipy import signal
from src.motion_analysis.MyVideoHelper2 import MyVideoHelper2
from abc import ABC, abstractmethod


class AbstractDenseOpticalFlow(MyVideoHelper2, ABC):
    """
    Represents a video of a sequence of images. It methods for estimating the optical flow
    are implemented with the Horn-Shunck method
    """
    def __init__(self, images):
        """
        Creates a MyVideo3 object from a sequence of two or more images of
        the same size
        
        params:
          -images: a list of two or more MyImage objects of the same size
        """
        
        MyVideoHelper2.__init__(self, images)    
    
    
    @abstractmethod
    def get_optical_flow_ROI(self, min_corner, max_corner, t, initial_displacements = None, max_iterations = 10000,
                             smoothness = 0.1, input_mask=None):
        """
        Given the upper-left corner and lower-right corner of 
        a region of interest and a time t, it calculates the
        optical flow
        
        params:
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] <= max_corner[0] < self.width
                * min_corner[1] <= max_corner[1] < self.height
            
            -t: it must satisfy 0 <= t < self.duration-1 
            
            -initial_displacements: a numpy array of size W x H x 2, where W and H are the width and height of the
            region of interest, used to initialize the algorithm. If not specefied, the algorithm starts with all zeros
            
            -max_iterations: maximum number of iterations that the algorithm can perform
            
            -smoothness: positive flow used to control the smoothness of the optical flow. The greater the smoother
        
        returns:
            a numpy array of floats of size W x H x 2, where W and H are the width and height of the
            region of interest, in which the values at (x, y, 0) and (x, y, 1) represent the x-component
            and y-component of the optical flow at pixel (x, y) relative to the ROI
        """
        pass
        
    
    