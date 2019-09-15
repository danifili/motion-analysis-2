import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

class BezierCurve:
  def __init__(self, start, end, control):
    self.start = np.array(start)
    self.end = np.array(end)
    self.control = np.array(control)
  
  def split_points(self, delta=1):
    length = self.estimate_length()
    num_steps = int(self.estimate_length() / delta)
    return np.array([self.get_point(self.find_t(length * step/num_steps)) for step in range(num_steps)])
  
  def estimate_length(self, num_steps=100, start_t=0, end_t=1):
    points = [self.get_point(start_t + (end_t - start_t)*step/num_steps) for step in range(num_steps+1)]
    distance = 0
    for idx in range(len(points)-1):
      distance += np.linalg.norm(points[idx] - points[idx+1])
    return distance
    
  def get_point(self, t):
    return (1-t)**2 * self.start + 2*(1-t)*t*self.control + t**2* self.end
  
  def find_t(self, sub_length, start_t=0, end_t=1, eps = 1e-6):
    mid_t = (start_t + end_t) / 2
    length_to_mid_t = self.estimate_length(start_t=start_t, end_t=mid_t)
    offset = sub_length - length_to_mid_t
    if np.abs(offset) < eps:
      return mid_t
    elif offset < 0:
      return self.find_t(sub_length, start_t=start_t, end_t=mid_t, eps=eps)
    else:
      return self.find_t(sub_length-length_to_mid_t, start_t=mid_t, end_t=end_t, eps=eps)
      

def read_csv_file(filename):
  with open(filename, 'r') as f:
    reader = csv.reader(f)
    array = []
    for row in reader:
      array.append(list(map(float,row)))
    return np.array(array)

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
    start = start_point * np.array([width/100, height/100]),
    control = control_point * np.array([width/100, height/100]),
    end = end_point * np.array([width/100, height/100])
  )

  points = bezier_curve.split_points(delta=1)

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

    


  



  







