from MyVideoHelper2 import MyVideoHelper2
import numpy as np
import functools
from timeit import default_timer as timer

from scipy import ndimage

class TimonerFreeman(MyVideoHelper2):
    """
    Represents a video of a sequence of images. Its methods for estimating the displacements
    of regions of interest are implemented with the method found in the Timoner-Freeman paper.
    """

    def __init__(self, images, Ex=None, Ey=None, Et=None):
        """
        params:
            -images: a numpy array of MyImage objects with the same with and height.
            
            -Ex: optional parameter. It is a 3-dimensional numpy array of shape W x H x 8, where
                W is the width of the images, H is the height of the images. The value at (x, y, t)
                represents the gradient with respect to the x-direction of the pixel (x, y) at time t.
                If not given, it will be computed as in the Timoner-Freeman paper.

            -Ey: optional parameter. It is a 3-dimensional numpy array of shape W x H x 8, where
                W is the width of the images, H is the height of the images. The value at (x, y, t)
                represents the gradient with respect to the y-direction of the pixel (x, y) at time t.
                If not given, it will be computed as in the Timoner-Freeman paper.

            -Et: optional parameter. It is a 3-dimensional numpy array of shape W x H x 8, where
                W is the width of the images, H is the height of the images. The value at (x, y, t)
                represents the gradient with respect to the time of the pixel (x, y) at time t.
                If not given, it will be computed as in the Timoner-Freeman paper.

        """
        MyVideoHelper2.__init__(self, images)
        self.__images = np.array(images)
        
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
        
                
        
    def get_displacement_ROI(self, min_corner, max_corner, t):
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

        
        returns:
            A 2-length numpy arrat representing the average displacement of region of interest from time 
            t to time (t+1). The first element represents the x direction and the y element represents
            the y direction
        
        """
#        return np.array([-0, 0])
        
        return self._gradient_based_estimate(min_corner, max_corner, t)


    def _gradient_based_estimate(self, min_corner, max_corner, time, offset = (0, 0), correlation = 50):
        """
		Computes the gradient based estimate displacement of a region of interest.

		params:
		  -min_corner: an integer tuple of length 2 representing the upper-left corner 
		   of the region of interest
		   It must satisfy:
                * 0 <= min_corner[0] < self.width
                * 0 <= min_corner[1] < self.height

		  -max_corner: an integer tuple of length 2 representing the lower-left corner 
		   of rthe egion of interest
		   It must satisfy:
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height

          -time: an integer representing the time of the displacement we are calculting.
           it must satisfy 0 <= time < self.duration-1

          -offset: an integer tuple of length 2 representing the offset by which the
           image is being shifted

        return:
          A numpy array length 2 representing the gradient based estimate displacement
          of region of interest

        """
        
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        sx, sy = offset
        
        
        #get gradients with respect to x and y
        gx = self.__Ex[x_min:x_max+1, y_min:y_max+1, time].ravel()
        gy = self.__Ey[x_min:x_max+1, y_min:y_max+1, time].ravel()
        gt = self.__Et[x_min:x_max+1, y_min:y_max+1, time].ravel()
        
        A = np.array([[np.dot(gx, gx), np.dot(gx, gy)],
					 [np.dot(gx, gy), np.dot(gy, gy)]])

        b = np.array([np.dot(gx, gt), np.dot(gy, gt)])

        
        return -1 * self._solve(A, b, correlation)

    
    def get_eigenvalues(self, min_corner, max_corner, time):
        """
        Get the eigenvalues of the matrix of gradients used to estimate the displacement
        in the given region of interest.
        
        params:
		  -min_corner: an integer tuple of length 2 representing the upper-left corner 
		   of the region of interest
		   It must satisfy:
                * 0 <= min_corner[0] - offset[0] < self.width
                * 0 <= min_corner[1] - offset[1] < self.height

		  -max_corner: an integer tuple of length 2 representing the lower-left corner 
		   of rthe egion of interest
		   It must satisfy:
                * min_corner[0] < max_corner[0] < self.width
                * min_corner[1] < max_corner[1] < self.height
            
          -time: an integer; it must satisfy 0 <= time < self.duration-1
         
        returns:
            A numpy array of floats of length 2 containing the eigenvalues of the matrix
            used to compute the gradient-based displacements
            
        """
        x_min, y_min = min_corner
        x_max, y_max = max_corner
        
        #get gradients with respect to x and y
        gx = self.__Ex[x_min:x_max+1, y_min:y_max+1, time].ravel()
        gy = self.__Ey[x_min:x_max+1, y_min:y_max+1, time].ravel()
        

        A = np.array([[np.dot(gx, gx), np.dot(gx, gy)],
					 [np.dot(gx, gy), np.dot(gy, gy)]])

        u, s, v = np.linalg.svd(A)
        return s


    def _solve(self, A, b, correlation):
        """
        Find the solution of the equation Ax = b by using SVD
        
        params:
            -A: 2-dimensional numpy array of size NxN
            -b: numpy array of lenght Nx1
            -correlation: a positive float which will ignore singular values 
                          which ratio to the maximum singular value is bigger
        
        return:
            A numpy array of size N with the solution. If matrix is singular,
            the zero vector will be returned
            
        """
        try:
            
            U, s, V = np.linalg.svd(A, full_matrices = True)

            s_inv = np.diag([1/x if x * correlation > max(s) else 0 for x in s])
            A_pseudoinverse = np.dot(np.transpose(V), np.dot(s_inv, np.transpose(U)))
            return np.dot(A_pseudoinverse, b)

        except np.linalg.linalg.LinAlgError:
            return np.zeros(A.shape[0])



