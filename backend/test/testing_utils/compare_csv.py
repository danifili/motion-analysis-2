import sys
import numpy as np
from src.utils.file_reader import read_csv_file
from src.motion_analysis.Plot import Plot
from src.motion_analysis.MyImage import MyImage
import matplotlib.pyplot as plt

if __name__ == "__main__":
  f1 = sys.argv[1]
  f2 = sys.argv[2]
  f3 = sys.argv[3]

  matrix1 = read_csv_file(f1)
  matrix2 = read_csv_file(f2)
  background_image = MyImage(f3)

  difference = np.abs((matrix1 - matrix2))
  plt.imshow(background_image[0:background_image.width, 0:background_image.height].T, cmap='gray',vmin=0,vmax=255)
  plt.imshow(difference.T, cmap='hot', interpolation='nearest', alpha=1)
  plt.colorbar()
  plt.show()