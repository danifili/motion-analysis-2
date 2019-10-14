from src.motion_analysis.MyImage import MyImage
from src.motion_analysis.MyVideoFinal import MyVideo
from src.motion_analysis.Plot import Plot
from src.utils.compare_images import plot_matrices_differences
from src.utils.compress_matrix import compress_matrix
from src.utils.file_reader import read_csv_file

import numpy as np
import matplotlib.pyplot as plt


def reconstruct_A_and_B(amps, phases):
  assert amps.shape == phases.shape
  As = amps * np.cos(phases)
  Bs = As * np.tan(phases)
  return As, Bs

def compress(amps, phases):
  As, Bs = reconstruct_A_and_B(amps, phases)
  compressed_As, compressed_Bs = compress_matrix(As, scale=4), compress_matrix(Bs, scale=4)

  new_amps = np.sqrt(compressed_As * compressed_As + compressed_Bs * compressed_Bs)
  new_phases = np.arctan2(compressed_Bs, compressed_As)
  new_phases[new_phases < 0] += 2 * np.pi

  return new_amps / 4, new_phases % (2 * np.pi)


if __name__ == "__main__":
  compressed_image_root = "./backend/testing_samples/"
  full_image_root = "./backend/testing_samples/full_"

  suffixes = ["amplitudes_x.csv", "amplitudes_y.csv", "phases_x.csv", "phases_y.csv"]

  compressed_image_results_files = [compressed_image_root + suffix for suffix in suffixes]
  full_image_results_files = [full_image_root + suffix for suffix in suffixes]

  compressed_image_results = [read_csv_file(filepath) for filepath in compressed_image_results_files]
  print ("read compressed images results")
  full_image_results = [read_csv_file(filepath) for filepath in full_image_results_files]
  print ("read full images results")

  amps_x, phases_x = full_image_results[0], full_image_results[2]
  amps_y, phases_y = full_image_results[1], full_image_results[3]

  compressed_amps_x, compressed_phases_x = compress(amps_x, phases_x)
  print ("compressed full image amps and phases in the x direction")
  compressed_amps_y, compressed_phases_y = compress(amps_y, phases_y)
  print ("compressed full image amps and phases in the y direction")

  post_compression_results = [compressed_amps_x, compressed_amps_y, compressed_phases_x, compressed_phases_y]

  background_image = MyImage("./backend/testing_samples/test_image.png", compress_image=True)

  for i in range(4):
    plt.figure()
    plt.title(suffixes[i])
    plot_matrices_differences(compressed_image_results[i], post_compression_results[i], background_image, mode = "percent" if i < 2 else "phase")
  plt.show()


