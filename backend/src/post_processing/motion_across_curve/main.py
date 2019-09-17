import numpy as np
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

from src.post_processing.motion_across_curve.bezier_curve import BezierCurve
from src.utils.file_reader import read_csv_file

def get_new_amp_and_phase(long_vec, amp_vec, phase_vec):
  theta = np.arctan2(long_vec[1], long_vec[0])
  Ax, Bx = amp_vec[0] * np.cos(phase_vec[0]), amp_vec[0] * np.sin(phase_vec[0])
  Ay, By = amp_vec[1] * np.cos(phase_vec[1]), amp_vec[1] * np.sin(phase_vec[1])
  
  new_Ax = Ax * np.cos(theta) + Ay * np.sin(theta)
  new_Bx = Bx * np.cos(theta) + By * np.sin(theta)
  new_Ay = -Ax * np.sin(theta) + Ay * np.cos(theta)
  new_By = -Bx * np.sin(theta) + By * np.cos(theta)

  new_amp_x = np.sqrt(new_Ax ** 2 + new_Bx ** 2)
  new_amp_y = np.sqrt(new_Ay ** 2 + new_By ** 2)

  new_phase_x = np.arctan2(new_Bx, new_Ax) % (2 * np.pi)
  new_phase_y = np.arctan2(new_By, new_Ay) % (2 * np.pi)

  return np.array([new_amp_x, new_amp_y]), np.array([new_phase_x, new_phase_y])

if __name__ == "__main__":
  start_point = np.array(list(map(float, sys.argv[1:3])))
  end_point = np.array(list(map(float, sys.argv[3:5])))
  control_point = np.array(list(map(float, sys.argv[5:7])))

  files = [sys.argv[i] for i in range(7,11)]
  files.sort()

  ampX = read_csv_file(files[0])
  ampY = read_csv_file(files[1])
  phaseX = read_csv_file(files[2])
  phaseY = read_csv_file(files[3])

  print ("read csv files")

  assert ampX.shape == ampY.shape == phaseX.shape == phaseY.shape

  width, height = ampX.shape

  bezier_curve = BezierCurve(
    start_point = start_point * np.array([width/100, height/100]),
    control_point = control_point * np.array([width/100, height/100]),
    end_point = end_point * np.array([width/100, height/100])
  )

  points = bezier_curve.split_curve_equally(curve_segment_length=1)

  print("divided bezier curve")

  long_amps = []
  long_phases = []
  radial_amps = []
  radial_phases = []

  for i in tqdm(range(len(points)-1)):
    point = points[i]
    next_point = points[i+1]
    discrete_x, discrete_y = int(round(point[0])), int(round(point[1]))
    
    (long_amp, radial_amp), (long_phase, radial_phase) = get_new_amp_and_phase(
      long_vec = next_point - point,
      amp_vec = np.array([ampX[discrete_x, discrete_y], ampY[discrete_x, discrete_y]]),
      phase_vec = np.array([phaseX[discrete_x, discrete_y], phaseY[discrete_x, discrete_y]])
    )

    long_amps.append(long_amp)
    long_phases.append(long_phase)
    radial_amps.append(radial_amp)
    radial_phases.append(radial_phase)
  
  plt.figure("Long Amp")
  plt.plot(long_amps)
  plt.figure("Long Phase")
  plt.plot(long_phases)
  plt.figure("Radial Amp")
  plt.plot(radial_amps)
  plt.figure("Radial Phase")
  plt.plot(radial_phases)

  plt.show()

    


  



  







