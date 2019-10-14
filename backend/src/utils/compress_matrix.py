import numpy as np

def compress_matrix(matrix, scale):
  """
  Compress matrix by using DFT by a factor of scala in both x and y axis.
  
  matrix: two dimensional numpy array
  scale: int

  Returns: compressed matrix
  """
  width, height = matrix.shape
  fft_image = np.fft.fft2(matrix)
  fft_image[int(width/(2 * scale)):-int(width/(2 * scale)),int(height/(2 * scale)) :-int(height/(2 * scale))] = 0
  return np.fft.ifft2(fft_image).real[::scale, ::scale]