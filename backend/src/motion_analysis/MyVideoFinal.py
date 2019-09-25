#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:31:12 2017

@author: danifili
"""
from src.motion_analysis.MyVideoHelper2 import MyVideoHelper2
from src.motion_analysis.MyVideoOpticalFlow import HornShunck
from src.motion_analysis.TimonerFreeman import TimonerFreeman
from multiprocessing import Pool, Process
import time
from tqdm import tqdm
from numba import njit
import numpy as np
from tqdm import tqdm

@njit
def sinusoidal_fit(cumulative_displacements):
    """
    Given a region of interest that have sinusoidal movement with an expected frequency, it calculates the
    amplitude and phase for every pixel
    
    params:
        -min_corner: A tuple of size 2 representing the (x, y) coordinates
            of the upper-left corner of our region of interest.
            It must satisfy:
            * 0 <= min_corner[0] < self.width
            * 0 <= min_corner[1] < self.height
            
        -max_corner: A tuple of size 2 representing the (x, y) coordinates 
        of the lower-right corner of our region of interest.
        It must satisfy:
            * min_corner[0] < max_corner[0] < self.width
            * min_corner[1] < max_corner[1] < self.height
        
        -win_min: The minimum size of the regions used to calculate the initial displacements of this algorithm
        
        -win_max: The maximum size of the regions used to calculate the initial displacements of this algorithm
        
        -max_iterations: maximum number of iterations that the algorithm can perform
        
        -smoothness: positive flow used to control the smoothness of the optical flow. The closer to 0 the smoother
        
        -min_eigen: a positive float; the min acceptable ratio between the eigenvalues and the area
        of the region of interest of the matrix used to compute its displacement.
        
    returns: 
        a numpy array of size W x H x 4, where W and H are the width and height of the region of interest, in which
        the (x, y) entry is a numpy array of size 4 containing amplitude_x, amplitude_y, phase_x, phase_y
    """
    period = 8
    frequency = 2 * np.pi / period
    t = np.array([frequency * p for p in range(period)])
    sin_t = np.sin(t)
    cos_t = np.cos(t)
    cos_t_minus_1 = cos_t - 1
    M = np.array([list(sin_t), list(cos_t_minus_1)])
    P = np.dot(np.linalg.inv(np.dot(M, M.T)), M)

    width = cumulative_displacements.shape[1]
    height = cumulative_displacements.shape[2]

    data = np.zeros((width, height, 4), dtype=np.float64)
    single_cumulative_displacements_x = np.zeros(8, dtype=np.float64)
    single_cumulative_displacements_y = np.zeros(8, dtype=np.float64)

    for x in range(width):
        for y in range(height):
            for f in range(period):
                single_cumulative_displacements_x[f] = cumulative_displacements[f, x, y, 0]
                single_cumulative_displacements_y[f] = cumulative_displacements[f, x, y, 1]
            amp_x, phase_x = sinusoidal_fit_helper(single_cumulative_displacements_x, P)
            amp_y, phase_y = sinusoidal_fit_helper(single_cumulative_displacements_y, P)

            data[x, y, 0] = amp_x
            data[x, y, 1] = amp_y
            data[x, y, 2] = phase_x
            data[x, y, 3] = phase_y

    return data

@njit
def sinusoidal_fit_helper(cumulative_displacements, P):
    """
    Given the displacements in a certain direction for a pixel that have sinusoidal movement with a given frequency,
    it calculates the amplitude and phase of the sine wave
    
    params:
        -displacements: a numpy array containing the displacements in a certain direction of two consecutive images
        
    returns:
        a tuple in which the first entry is the amplitude of the sine wave and the second entry is the phase of
        the sine wave
    """

    projection = P.dot(cumulative_displacements)
    A = projection[0]
    B = projection[1]

    amp = np.sqrt(A * A + B * B)
    phase = np.arctan2(B, A)

    if phase < 0:
        phase += 2 * np.pi

    return amp, phase % (2 * np.pi)

class MyVideo(MyVideoHelper2):
    
    def __init__(self, images):
        MyVideoHelper2.__init__(self, images)
        
        #precompute matrices of gradients
        min_corner = (0, 0)
        max_corner = (self.width-1, self.height-1)
        
        image_shape = (self.width, self.height, self.duration-1)
        self.__Ex = np.zeros(image_shape)
        self.__Ey = np.zeros(image_shape)
        self.__Et = np.zeros(image_shape)

        print ("Computing gradients...")
        for t in tqdm(range(self.duration-1)):
            self.__Ex[:, :, t] = self._Ex(min_corner, max_corner, t)
            self.__Ey[:, :, t] = self._Ey(min_corner, max_corner, t)
            self.__Et[:, :, t] = self._Et(min_corner, max_corner, t)

        self.__video_displacement_ROI = TimonerFreeman(images, self.__Ex, self.__Ey, self.__Et)
        self.__video_optical_flow = HornShunck(images, self.__Ex, self.__Ey, self.__Et)
    
    def get_displacement_ROI(self, min_corner, max_corner, t, max_margin = 0):
        """
        Given the upper-left corner and lower-right corner of 
        a region of interest and a time t, it calculates the displacement 
        between time t and t+1 in such ROI
        
        params:
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: it must satisfy 0 <= t < self.duration-1
            
            -max_margin: a positive integer; the maximum margin used to compute the displacement
            in the region of interest.

        
        returns:
            A tuple of length 3 containing:
             -A 2-length numpy array representing the average displacement of region of interest from time 
              t to time (t+1). The first element represents the x direction and the y element represents
              the y direction

             -the two eingenvalues that come from the matrix used to compute the displacemets divide by
              the number of pixels used to compute this matrix
        """
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        roi_width = x_max - x_min + 1
        roi_height = y_max - y_min + 1
        
        #find the best margin
        x_margin = min(x_min, self.width-1-x_max)
        y_margin = min(y_min, self.height-1-y_max)
        margin = min(max_margin, x_margin, y_margin)
        
        new_min_corner = (x_min - margin, y_min - margin)
        new_max_corner = (x_max + margin, y_max + margin)
        
        #get the minimum eigenvalue of the region of interest to determine if there is enough contrast
        new_area = (2*margin + roi_width) * (2*margin + roi_height)
        displacements, (l1, l2) = self.__video_displacement_ROI.get_displacement_ROI(new_min_corner, new_max_corner, t)
        return displacements, l1 / new_area, l2 /new_area 

          
    def get_optical_flow_ROI(self, min_corner, max_corner, t, win_min = 5, win_max = 25, quality_level = 0.07, max_iterations = 100,  smoothness = 100):
        """
        Given the upper-left corner and lower-right corner of 
        a region of interest and a time t, it calculates the
        optical flow
        
        params:
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 < min_corner[0] < self.width
                * 0 < min_corner[1] < self.height
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
            -t: it must satisfy 0 <= t < self.duration-1 
            
            -win_min: The minimum size of the regions used to calculate the initial displacements of this algorithm
            
            -win_max: The maximum size of the regions used to calculate the initial displacements of this algorithm
            
            -max_iterations: maximum number of iterations that the algorithm can perform
            
            -smoothness: positive flow used to control the smoothness of the optical flow. The closer to 0 the smoother
            
            -min_eigen: a positive float; the min acceptable ratio between the eigenvalues and the area
            of the region of interest of the matrix used to compute its displacement.
        
        returns:
             a numpy array of floats of size W x H x 2, where W and H are the width and height of the
             region of interest, in which the values at (x, y, 0) and (x, y, 1) represent the x-component
             and y-component of the optical flow at pixel (x, y) relative to the ROI
        """
        print ("Computing displacements from frame " + str(t) + " to frame " + str(t+1))
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        roi_width = x_max - x_min + 1
        roi_height = y_max - y_min + 1
        
        #initialize initial displacements
        displacements_mask = np.zeros((roi_width, roi_height, 2))
        
        #get maximum margin
        max_margin = (win_max - win_min) // 2
        
        points = []
        displacements = []
        mask = np.zeros((roi_width, roi_height, 2))

        for x in tqdm(range(roi_width // win_min)):
            for y in range(roi_height // win_min):
                new_min_corner = np.array([win_min * x + x_min, win_min * y + y_min])
                new_max_corner = np.array([win_min * (x+1) + x_min-1, win_min*(y+1)+y_min-1])
                
                x_min_rel, y_min_rel = new_min_corner - np.array([x_min, y_min])
                x_max_rel, y_max_rel = new_max_corner - np.array([x_min, y_min])
                
                
                #get average displacements 
                avg_displacement, l1, l2 = self.get_displacement_ROI(new_min_corner, new_max_corner, t,\
                                                            max_margin)

                mask[(x_min_rel + x_max_rel)//2, (y_min_rel + y_max_rel)//2, :] = np.array([l1, l2])
                displacements_mask[(x_min_rel + x_max_rel)//2, (y_min_rel + y_max_rel)//2, :] = np.array(avg_displacement)

        mask[:, :, 0] = mask[:, :, 0]  / np.max(mask[:, :, 0])
        mask[:, :, 1] = mask[:, :, 1]  / np.max(mask[:, :, 1])

        for x in range((win_min-1)//2, roi_width, win_min):
            for y in range((win_min-1)//2, roi_height, win_min):

                if min(mask[x,y]) >= quality_level:
                    displacements.append(displacements_mask[x,y])
                    points.append([x,y])

            
        try:
            displacements = self._interpolate(roi_width, roi_height, np.array(points), np.array(displacements))
        except:
            raise ValueError("Value of input quality error " + str(quality_level) + " is too large. Try an smaller value")

        #use algorithm for calculating optical flow given the average displacements as the initial ones
        return self.__video_optical_flow.get_optical_flow_ROI(min_corner, max_corner, t, \
                initial_displacements = displacements, max_iterations = max_iterations, smoothness = smoothness, input_mask=mask)
    
    def get_cumulative_displacements(self, min_corner, max_corner, win_min = 5, win_max = 25, quality_level = 0.07, max_iterations = 100,  smoothness = 100):
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        roi_width = x_max - x_min + 1
        roi_height = y_max - y_min + 1

        directions = 2

        cumulative_displacements = np.zeros((self.duration, roi_width, roi_height, directions), dtype=np.float64)

        for t in range(self.duration-1):
            new_displacements = self.get_optical_flow_ROI(min_corner, max_corner, t, win_min=win_min, win_max=win_max, quality_level=quality_level,\
                                                    max_iterations=max_iterations, smoothness=smoothness)

            cumulative_displacements[t+1] = cumulative_displacements[t] + new_displacements

        return cumulative_displacements

if __name__ == "__main__":
    for i in tqdm(range(10)):
        sinusoidal_fit(np.random.rand(8, 2500, 2500, 2))
    