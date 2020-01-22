from PIL import Image as Img
import numpy as np
from scipy.interpolate import griddata

import mmap
import contextlib
import re
import struct
import datetime

from src.utils.compress_matrix import compress_matrix

class MyImage(object):
    """
    Represents a matrix of the pixels of an image
    """
    FIXED_PATTERN_NOISE = "fixed pattern noise"
    SHOT_NOISE = "shot noise"
    NO_NOISE = "no noise"
    LINEAR = "linear"
    FFT = "fft"
    _MAX_COLOR = 255
    
    
    """
    Abstraction function:
        AF(pixels): An image in which the the brightness in the pixel (x, y) is pixels[x, y].
        The x-coordinates increase to right and the y-coordinates increase
        downwards.
    
    Representation Invariant:
        pixels is a 2-dimensional numpy array. pixels[x, y] is a integer from 0 to 255
    
    Safety from rep exposure:
        Field is private and never mutated nor reassigned. When created, a defensive copy of the
        array is made
    
    """
    
    def __init__(self, file_path, compress_image=False):
        """
        Creates an MyImage object from the file path of the image we want to represent
        params:
            -file_path: string of the directory of the image we want to represent
        """
        
        
        

        #f = open(file_path,'r')
        
        
        #print(self.__pixels[0][1])

        with open(file_path, 'r') as f:
            f.seek(0,2) #position at end of file
            filesize = f.tell()
            f.seek(0,0) # back to beginning
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            if m.read(2) != b'nD':        
                pic = Img.open(file_path).convert('L')
                self.__pixels = np.transpose(np.array(pic))
            else:                
                pattern = re.compile(b'nD (\d+) (\d+)\n# 2D (\w+) (\d+) (\d+)') #http://docs.python.org/2/library/re.html
                (dataoffset, notesoffset, pixeltype, width, height) = [t(s) for t,s in zip((int, int, bytes, int, int),re.search(pattern,m).groups())]
                print (file_path, ' filesize: ', filesize, ' dataoffset: ', dataoffset, ' notesoffset: ', notesoffset, ' pixeltype: ', pixeltype, ' width: ', width, ' height: ', height)
                if dataoffset > notesoffset | notesoffset > filesize:
                    print("Bad nD file: file too short: ", file_path)
                    #exit with error
                if (pixeltype) == b'F32':
                    if notesoffset - dataoffset != 4*width*height:
                        print("2D file has wrong length: ", file_path)
                        #exit with error
                else:
                    print("Unhandled nD type: ", pixeltype.decode()) #Could add support for other types ...
                    #exit with error

                    
                self.__pixels = np.zeros((width, height),dtype=np.float32) #initialize array of F32 values
                m.seek(dataoffset,0)
                
                for j in range(0,height):
                    for i in range(0,width):                            
                        self.__pixels[i][j] = struct.unpack('<f',m.read(4))[0] 

                

                if True:  #choose to dscale or not
                    ############Scale from F32 for rest of program                
                    firstpointzeroflag = False 
                    if self.__pixels[0][0] == 0: #In my datasets the first [0,0] pixel is always 0 for some reason. Not sure why, but can exclude.
                        npmin = np.min(self.__pixels[1:]) #exclude first point to increase dynamic range and reduce quantization in conversion
                        npmax = np.max(self.__pixels[1:])
                        firstpointzeroflag = True
                    else:
                        npmin = np.min(self.__pixels)
                        npmax = np.max(self.__pixels)
                    #print("Downsampling from F32 (", npmin, "-", npmax, ") to Int8 (0-255)...")
                    print("Scaling F32 (" + str(npmin) + "-" + str(npmax) + ") to F32 (0.0-" + str(self._MAX_COLOR) + ".0)...")
                    ratio = self._MAX_COLOR/(npmax - npmin) #this is wrong as it should look at min and max for whole set
                    self.__pixels = (self.__pixels-npmin) *ratio
                    
                    if firstpointzeroflag == True:
                        self.__pixels[0][0] = 0 #first point was 0, so put it back if it changed

                    #commented 2 lines below out as it seems like the program can handle floats between 0-255?
                    #So output is currently a 32-bit float scaled into 0-255 without the below line. However, the scaling changes between images right now.
                        
                    #self.__pixels = np.rint(self.__pixels).astype(np.int) #round and convert to int
                        
                    ############


                
                ################## Save the image as nD
                #test our python nD writing powers
                # writefile=file_path + ".new"
                # with open(writefile, 'bw') as w:                   
                #     prefix = b'# 2D F32 ' + str(width).encode() + b' ' + str(height).encode() + b'\n=date=' + datetime.datetime.now().strftime("%Y-%m-%d").encode()
                #     prefix += b'\n=time=' +  datetime.datetime.now().strftime("%H:%M:%S").encode() + b'\n=command=\n' #could add commandline w/ params here before \n
                    
                #     nprefix = len(prefix)
                #     dataoffset = nprefix
                #     datasize = 4*width*height
                #     notesoffset = dataoffset+datasize
                #     notes=b'This file was created with motion analysis'
                #     header = b'nD ' + str(dataoffset).encode() + b' ' + str(notesoffset).encode() + b'\n'
                #     while dataoffset != len(header) + nprefix:
                #         dataoffset = len(header) + nprefix
                #         notesoffset = dataoffset + datasize
                #         header = b'nD ' + str(dataoffset).encode() + b' ' + str(notesoffset).encode() + b'\n'

                    
                #     w.write(header + prefix)
                #     for j in range(0, height):
                #         for i in range(0, width):
                #             w.write(struct.pack('f', self.__pixels[i][j]))
                #     w.write(notes)
                #     print("Wrote nD image to ", writefile)
                ###################################
                
                                
            m.close()
        if compress_image:
            print ("compressing image")
            self.__pixels = compress_matrix(self.__pixels, scale=4)
    
    #define width property
    def _width(self):
        return self.__pixels.shape[0]
    
    width = property(_width)
    
    #define height property
    def _height(self):
        return self.__pixels.shape[1]
    
    height = property(_height)
    
    def __getitem__(self, coordinates):
        """
        Returns the brightness of pixels given their coordinates
        
        params:
            -coordinates: A two dimensional tuple (x, y) of integers in which:
                        * 0 <= x < self.width and x = 0 represents left of image
                        * 0 <= y < self.height and y = 0 represents the top of image
                
                        or a tuple of slices (x0:x1, y0: y1) in which:
                            * 0 <= x0 < x1 <= self.width
                            * 0 <= y0 < y1 <= self.height
        
        returns:
            If coordinates is a tuple (x, y), an integer between 0 and 255 (inclusive) 
            which represents the brightness of pixel will be return.
            
            Otherwise, return a numpy array of shape (x1-x0, y1-y0) in which the (x, y)
            component corresponds to an integer from 0 to 255 (inclusive) that represents the
            brightness at pixel (x0+x, y0+y) in the image
                
        """
        #fail fast if preconditions are not met
        assert type(coordinates) is tuple
        assert len(coordinates) == 2
        assert type(coordinates[0]) is type(coordinates[1])
        
        x, y = coordinates
        
        if type(x) is int:
            #put x and y in the correct range
            x, y = max(0, x), max(0, y)
            x, y = min(x, self.width-1), min(y, self.width-1)
            return self.__pixels[x,y]
        
        
        if type(x) is slice:
            x0, x1, y0, y1 = x.start, x.stop, y.start, y.stop
            assert x0 is not None and x1 is not None and y0 is not None and y1 is not None, "all 4 points in slice must be specified"
            pixels = np.zeros((x1-x0, y1-y0))
            
            new_x = slice(max(0, x0), min(self.width, x1))
            new_y = slice(max(0, y0), min(self.height, y1))
            
            new_x_rel = slice(max(0, -x0), x1-x0 - max(0, x1-self.width))
            new_y_rel = slice(max(0, -y0), y1-y0 - max(0, y1-self.height))
            
            pixels[new_x_rel, new_y_rel] = self.__pixels[new_x, new_y]
            
            return pixels
            
        
        raise TypeError("coordinates must be a tuple of length 2 of integers or slices")

    @staticmethod
    def image_from_matrix(pixels, file_path, noise = NO_NOISE):
        """
        Returns a MyImage object given a numpy matrix of integers and creates and stores its 
        equivalent BMP image on the given file path. 
        Integer values of matrix must be between 0 and 255 (inclusive). 
        
        params:
            pixels: numpy matrix of floats between 0 and 255 (inclusive)
            file_path: string representing name of file and directory. It must finish in .bmp
            
            noise: an MyImage flag. It can be either MyImage.NO_NOISE, MyImage.SHOT_NOISE or MyImage.FIXED_PATTERN_NOISE.
        
        returns:
            MyImage object stored in file_path
        """
        
        width, height = pixels.shape
        im = Img.new("L", (width, height))
        picture = im.load()
        
        for x in range(width):
            for y in range(height):
                new_pixel = int(round(pixels[x,y]))
                if noise:
                    new_pixel = MyImage._add_noise(new_pixel, noise)
                picture[x,y] = new_pixel
        
        im.save(file_path, "BMP")
        im.close()


        ################## Save the image as nD
        writefile=file_path + ".nD"
        with open(writefile, 'bw') as w:                   
            prefix = b'# 2D F32 ' + str(width).encode() + b' ' + str(height).encode() + b'\n=date=' + datetime.datetime.now().strftime("%Y-%m-%d").encode()
            prefix += b'\n=time=' +  datetime.datetime.now().strftime("%H:%M:%S").encode() + b'\n=command=\n' #could add commandline w/ params here before \n
            
            nprefix = len(prefix)
            dataoffset = nprefix
            datasize = 4*width*height
            notesoffset = dataoffset+datasize
            notes=b'This file was created with motion analysis'
            header = b'nD ' + str(dataoffset).encode() + b' ' + str(notesoffset).encode() + b'\n'
            while dataoffset != len(header) + nprefix:
                dataoffset = len(header) + nprefix
                notesoffset = dataoffset + datasize
                header = b'nD ' + str(dataoffset).encode() + b' ' + str(notesoffset).encode() + b'\n'

            
            w.write(header + prefix)
            for j in range(0, height):
                for i in range(0, width):
                    w.write(struct.pack('f', pixels[i][j]))
            w.write(notes)
            print("Wrote nD image to ", writefile)
        ###################################
        
        
        return MyImage(file_path)

    @staticmethod
    def image_from_function(f, width, height, file_path, noise = NO_NOISE):
        """
        Returns a MyImage object given a function f, the width and height of
        the image, and creates and stores its equivalent BMP image in the given file_path.
        
        params:
            f: a function that maps the integer coordinates (x, y) for
            0 <= x < width and 0 <= y < height to an float between 0 and 255
            
            width: a positive integer representing the width of the image
            
            height: a positive integer representing the height of the image
            
            file_path: string representing name of file and directory. 
            It must finish in .bmp
            
            noise: a MyImage flag. It can be either MyImage.NO_NOISE, MyImage.SHOT_NOISE or MyImage.FIXED_PATTERN_NOISE.

            
        returns:
            MyImage object of width and height given by the inputs such that
            f(x, y) rounded to the closest integer
            is the brightness of the pixel of coordinate (x, y). This
            image will be stored in the file_path given in the input.
        """
        
        assert width > 0
        assert height > 0
        
        #initialize matrix of pixels
        pixels = np.zeros((width, height))
        
        for x in range(width):
            for y in range(height):
                #get value of pixel
                pixel = int(round(f(x, y)))
                
                #assert that pixel is the correct range
                assert 0 <= pixel <= MyImage._MAX_COLOR
                
                #update brightness in pixel (x,y)
                pixels[x,y] = pixel
        
        return MyImage.image_from_matrix(pixels, file_path, noise)
    

    def shift_image(self, u, v, file_path, mode = LINEAR):
        """
        Create a shifted version of the current image given the displacement functions 
        
        params:
            u: a function whose inputs are the floats x, y representing x-coordinate, y-coordinate
            and returns a float representing the velocity in the x-direction
            
            v: a function whose inputs are the floats x, y representing x-coordinate, y-coordinate
            and returns a float representing the velocity in the y-direction
            
            file_path: string representing name of file and directory. 
            It must finish in .bmp.

            mode: a MyImage flag. It can be either MyImage.LINEAR, which will use a linear interpolator or
            MyImage.FFT, which will use FFT. If MyImage.FFT is used, then v(x, y1) = v(x, y2) and 
            u(x, y1) = u(x, y2) for every x, y1 and y2.
            
        returns: a MyImage object resulting from applying the given motion to this image
        """
        if mode == MyImage.LINEAR:
            return self._linear_shift_image(u,v,file_path)
        
        if mode == MyImage.FFT:
            return self._fft_shift_image(u,v,file_path)
        
        raise TypeError("invalid mode")

    def _linear_shift_image(self, u, v, file_path):
        shifted_pixels = np.zeros((self.width * self.height, 2))
        
        for x in range(self.width):
            for y in range(self.height):
                i = x * self.height + y
                dx = u(x, y)
                dy = v(x, y)
                shifted_pixels[i] = np.array([x+dx, y+dy])
        
        grid_x, grid_y = np.mgrid[0:self.width, 0:self.height]
        
        
        image = griddata(shifted_pixels, np.copy(self.__pixels).flatten(), (grid_x, grid_y), fill_value = 0)
        
        return MyImage.image_from_matrix(image, file_path)
    
    def _fft_shift_image(self, u, v, file_path):
        pixels = np.zeros((self.width, self.height))
        for x in range(self.width):
            meshy = np.zeros(self.height)
            for y in range(self.height):
                meshy[y] = (y - self.height/2) * v(x, y)
            
            dy = np.exp(-1j*2*np.pi*(meshy/self.height)) 
            
            column = np.copy(self.__pixels[x])
            column = np.fft.fftshift(np.fft.fft(column))
            column = np.multiply(column, dy)
            column = np.fft.ifft(np.fft.ifftshift(column))
            
            pixels[x] = np.copy(column.real)

        return MyImage.image_from_matrix(pixels, file_path)

    @staticmethod 
    def _add_noise(pixel, noise):
        """
        Adds noise to a pixel.
        
        params:
            pixel: the brightness of a certain pixel. It must be an integer from 0 to 255
            noise: a MyImage flag. It can be either MyImage.NO_NOISE, MyImage.SHOT_NOISE or MyImage.FIXED_PATTERN_NOISE.
        
        returns:
            a integer between 0 and 255 representing the new brightness after applying noise
        """
        if noise == MyImage.NO_NOISE:
            return pixel
        
        if noise == MyImage.SHOT_NOISE:
            return MyImage._add_shot_noise(pixel)
        
        if noise == MyImage.FIXED_PATTERN_NOISE:
            return MyImage._add_fixed_pattern_noise(pixel)
        
        raise TypeError("noise input is invalid")
    
    @staticmethod
    def _add_fixed_pattern_noise(pixel, variance=0.001):
        """
        Adds fixed pattern noise to a pixel
        
        params:
             pixel: the brightness of a certain pixel. It must be an integer from 0 to 255
             
             variance: the variance of the gaussian used to compute the noise
            
        returns: 
            an integer between 0 and 255 representing the new brightness after applying noise
        """
        mean = 1
       
        #update value of pixels
        noisy_pixel = pixel * np.random.normal(mean, variance)
        noisy_pixel = max(0, noisy_pixel)
        noisy_pixel = min(noisy_pixel, MyImage._MAX_COLOR)

            
        return int(round(noisy_pixel))
    
    @staticmethod
    def _add_shot_noise(pixel, scalar = 10):
        """
        Adds shot noise to a pixel
        
        params:
             pixel: the brightness of a certain pixel. It must be an integer from 0 to 255
             scalar: a positive integer used to scale the pixel brightness for estimating the
                     mean of the poison process. 
            
        returns: 
            an integer between 0 and 255 representing the new brightness after applying noise
        """
        mean = pixel * scalar
        noisy_pixel = np.random.poisson(mean) / scalar
        noisy_pixel = min(noisy_pixel, MyImage._MAX_COLOR)
        noisy_pixel = max(noisy_pixel, 0)
        
        return int(round(noisy_pixel))
    
    def show_image(self):
        pil_image = Img.fromarray(self.__pixels.T)
        pil_image.show()
            
    def get_pixels(self):
        #deep copy field for avoiding aliasing
        return np.array(self.__pixels)
    
if __name__ == "__main__":
    image = MyImage("backend/testing_samples/test_image.png")
    compressed_image = image.compress_image(4, "backend/testing_samples/test_compress_image.png")
    #compressed_image.show_image()

    print(compressed_image.width, compressed_image.height)


    #image = MyImage.image_from_matrix(image[300:1350, 500:1050], "tecta_art7_0.png")
    
    
#
#    k = 2 * np.pi / 500
#    w = 2 * np.pi / 8
#    
#    Ay = lambda x, y: np.exp((30 - x) /250 )
#    u = lambda x, y, t: 0
#    v = lambda x, y, t: Ay(x, y) * np.exp(complex(0, k*x - w*t)).real
#    
#    
#    for t in range(1, 8):
#        u_t = lambda x, y: u(x, y, t) - u(x, y, 0)
#        v_t = lambda x, y: v(x, y, t) - v(x, y, 0)
#        
#        image.shift_image(u_t, v_t, "tecta_art7_" + str(t) + ".png")
    
