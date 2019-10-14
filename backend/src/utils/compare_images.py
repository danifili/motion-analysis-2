import matplotlib.pyplot as plt
import numpy as np

def plot_matrices_differences(matrix1, matrix2, background_image, mode="unit"):
  difference = np.abs((matrix1 - matrix2))
  if mode == "unit":
    pass
  elif mode == "percent":
    matrix_abs_1 = np.abs(matrix1)
    matrix_abs_1[matrix_abs_1 == 0] = 1e-6
    difference = difference / matrix_abs_1 * 100
    difference[difference > 100] = 100 #clipping percentage values
  elif mode == "log":
    difference = np.log(difference)
  elif mode == "phase":
    difference[difference >= np.pi] = 2*np.pi - difference[difference >= np.pi]
  elif mode == "sqrt":
    difference = np.sqrt(difference)
  else:
    raise RuntimeError("wrong mode input value")
  plt.imshow(background_image[0:background_image.width, 0:background_image.height].T, cmap='gray',vmin=0,vmax=255)
  plt.imshow(difference.T, cmap = 'hot', interpolation='nearest', alpha=1)
  plt.colorbar()