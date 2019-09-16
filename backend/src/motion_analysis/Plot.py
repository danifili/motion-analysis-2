#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 13:53:13 2017

@author: danifili
"""
import matplotlib.pyplot as plt
# import cv2
import numpy as np
#import scipy.cluster.vq as kmeans
# from sklearn.cluster import KMeans


class Plot(object):
    
    @staticmethod
    def plot(video, displacements, min_corner, max_corner, time, k = 5, scale = 0.1, color = 'black'):
        """
        Plots the optical flow of a region of interest as a vector field
        
        params:
            -video: a MyVideo object
            
            -displacements: A 3-dimensional numpy array of floats containing the displacements
            of a certain region of interest, in which the value at (x, y, 0) represents the sdisplacement 
            in the x-direction and (x, y, 1) in the y-direction
            
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < self.width-2
                * 0 <= min_corner[1] < self.height-2
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < self.width-1
                * min_corner[1] < max_corner[1] < self.height-1
            
            -time: the the displacements are being calculated
            
            -scale: A float indicating the scale used to plot the vector field. For instance,
            if plot=True and scale=0.1, the arrows will be 10 times bigger.
            
            -k: A positive integer whose square is inversively proportional to the number of arrows plotted from
            the vector field.
            
            -color: a string; it indicates the color of the arrows of the vector field
        
        returns:
            nothing
        """
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        roi_width, roi_height = displacements.shape[:2]
        plt.imshow(video[x_min:x_max+1, y_min:y_max+1, time].T, cmap='gray',vmin=0,vmax=255)
        #get coordinates of arrows
        X, Y = [[k * x for x in range(roi_width//k)] for y in range(roi_height//k)], \
                 [[k * y for x in range(roi_width//k)] for y in range(roi_height//k)]
                 
        #get dimension of arrows 
        U, V = [[displacements[k * x, k * y, 0] for x in range(roi_width//k)] for y in range(roi_height//k)], \
                 [[- displacements[k * x, k * y, 1] for x in range(roi_width//k)] for y in range(roi_height//k)]
                  
        plt.quiver(X, Y, U, V, scale=scale, units="xy", color=color)
    
    @staticmethod
    def scatter_plot(image, min_corner, max_corner, x, y, color):
        """
        Given an image and a region of interest, a set of points and a color, it creates a scatter plot
        with the given points.

        params:
            -image: A MyImage object.
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < image.width-1
                * 0 <= min_corner[1] < image.height-1
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < image.width
                * min_corner[1] < max_corner[1] < image.height

            -x: A numpy array of integers containing the x-coordinates of the points to be plotted. x[i] corresponds to the
                x-coordinate of the i-th point, it is relative to the region of interest and it satisfies:
                    * min_corner[0] <= x[i] <= max_corner[0]

            -y: A numpy array of integers containing the y-coordinates of the points to be plotted. It must have the same length
                as x. y[i] corresponds to the y-coordinate of the i-th point, it is relative to the region of interest and it satisfies:
                    * min_corner[1] <= y[i] <= max_corner[1]

            -color: The color of the points to be plotted.

        returns:
            nothing
        """
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        plt.imshow(image[x_min:x_max+1, y_min:y_max+1].T)
        plt.scatter(x, y, c=color, s=1)
    
    
    @staticmethod
    def scalar_heat_map(video, min_corner, max_corner, time, scalars, alpha=1):
        plt.imshow(video[min_corner[0]: max_corner[0]+1, min_corner[1]:max_corner[1]+1, time].T, cmap='gray',vmin=0,vmax=255)
        plt.imshow(scalars.T, cmap='hot', interpolation='nearest', alpha=alpha)
        plt.colorbar()
    
    @staticmethod
    def vector_heat_map(video, min_corner, max_corner, time, flow, angle, alpha=1):
        """
        Given a matrix containing the optical flow of a region of interest, it plots
        a heat map of the magnitude of the flow in a given direction
        
        params:
            -video: a MyVideo object
            
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < self.width-2
                * 0 <= min_corner[1] < self.height-2
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < self.width-1
                * min_corner[1] < max_corner[1] < self.height-1
            
            -time: the the displacements are being calculated
            
            -flow:
                A 3-dimensional numpy array of shape W x H x 2, where the (x, y, 0) and (x, y, 1)
                entries correspond to the x and y components of the optical flow at pixel (x, y)
                relative to the region of interest
            
            -angle:
                A float representing direction. It must be in radians
        
        returns:
            None
        """
        plt.imshow(video[min_corner[0]: max_corner[0]+1, min_corner[1]:max_corner[1]+1, time].T, cmap='gray',vmin=0,vmax=255)
        A = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
        
        roi_width, roi_height, directions = flow.shape
        
        new_flow = np.array([[np.dot(A, flow[x, y]) for y in range(roi_height)] for x in range(roi_width)])
        
        plt.imshow(new_flow[:, :, 0].T, cmap='hot', interpolation='nearest', alpha=alpha)
        plt.colorbar()
    
    @staticmethod
    def phase_heat_map(video, min_corner, max_corner, time, phases, alpha=1):
        """
        Given a matrix containing phases in radians, it plots a hsv map
        
        params:
            -video: a MyVideo object
            
            -min_corner: A tuple of size 2 representing the (x, y) coordinates
             of the upper-left corner of our region of interest.
             It must satisfy:
                * 0 <= min_corner[0] < self.width-2
                * 0 <= min_corner[1] < self.height-2
                
            -max_corner: A tuple of size 2 representing the (x, y) coordinates 
            of the lower-right corner of our region of interest.
            It must satisfy:
                * min_corner[0] < max_corner[0] < self.width-1
                * min_corner[1] < max_corner[1] < self.height-1
            
            -time: the the displacements are being calculated
            
            -phases:
                A 2-dimensional numpy array of floats
        
        returns:
            None
        """
        plt.imshow(video[min_corner[0]: max_corner[0]+1, min_corner[1]:max_corner[1]+1, time].T, cmap='gray',vmin=0,vmax=255)
        plt.imshow(phases.T, cmap='hsv', interpolation='nearest', alpha=alpha, vmin=0, vmax=np.pi*2)
        plt.colorbar()
        

    @staticmethod
    def segmentation(image, flow, min_corner, max_corner, alpha, beta, k):
        colors = Plot.generate_random_colors(k)
        width, height = flow.shape[:2]
        vectors = np.zeros((width * height, 5))
        
        for x in range(width):
            for y in range(height):
                i = height * x + y
                dx, dy = flow[x, y]
                vectors[i] = np.array([x, y, alpha*dx, alpha*dy, beta*image[x, y]])
        
        kmeans = KMeans(k)
        labels = kmeans.fit_predict(vectors)

        
        clusters = [[[], []] for i in range(k)]
        
        for i in range(len(labels)):
            cluster = labels[i]
            x, y, adx, ady, bi = vectors[i]
            clusters[cluster][0].append(x)
            clusters[cluster][1].append(y)
        
        for i in range(k):
            x, y = clusters[i]
            color = colors[i]
            Plot.scatter_plot(image, min_corner, max_corner, x, y, color)
     
    @staticmethod
    def generate_random_colors(k):
        colors = np.zeros((k, 1, 3))
        for i in range(k):
            r = np.random.random()
            g = np.random.random()
            b = np.random.random()
            colors[i] = np.array([[r, g, b]])
        
        return colors
        
        
                
        
        
        
        
                
        
        
        

    
    