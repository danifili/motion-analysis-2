from src.utils.file_reader import read_csv_file, write_csv_file
from src.post_processing.motion_across_curve.bezier_curve import BezierCurve
from src.post_processing.motion_across_curve.curve_utils import *

import numpy as np
import matplotlib.pyplot as plt

class DisplacementThroughCurve:
  def __init__(self, start_point, end_point, control_point):
    self.start_point = start_point
    self.end_point = end_point
    self.control_point = control_point
  
  def _load_displacements_through_curve(self, displacementX_file, displacementY_file, cacheLoc = None):
    displacementsX = read_csv_file(displacementX_file)
    displacementsY = read_csv_file(displacementY_file)
    width, height = displacementsX.shape

    points = get_points_from_GUI(
      start_point=self.start_point,
      control_point=self.control_point,
      end_point=self.end_point,
      width=width,
      height=height,
      cacheLoc=cacheLoc
    )
    
    long_displacements_list = []
    radial_displacements_list = []
    for i in range(len(points)-1):
      current_point, next_point = points[i], points[i+1]
      discreteX, discreteY = clamp(int(round(current_point[0])), 0, width-1), clamp(int(round(current_point[1])), 0, height-1)
      long_dis, radial_dis = change_basis(
        reference_long_vec=next_point-current_point,
        vecX=displacementsX[discreteX, discreteY],
        vecY=displacementsY[discreteX, discreteY]
      )
      long_displacements_list.append(long_dis)
      radial_displacements_list.append(radial_dis)
    
    return np.array(long_displacements_list), np.array(radial_displacements_list)
  
  def plot_displacements(self, displacementX_file, displacementY_file, out_file_long, out_file_radial, cacheLoc = None):
    displacements_long, displacements_radial = \
      self._load_displacements_through_curve(
        displacementX_file=displacementX_file,
        displacementY_file=displacementY_file,
        cacheLoc=cacheLoc
      )

    plt.figure()

    plt.xlabel("Pixels")
    plt.ylabel("Displacement in Pixels")
    plt.plot(displacements_long)
    plt.hlines(0, 0, len(displacements_long)-1, linestyles="dashed")
    plt.ylim([-1.2,1.2])
    plt.savefig(out_file_long)

    plt.figure()

    plt.xlabel("Pixels")
    plt.ylabel("Displacement in Pixels")
    plt.plot(displacements_radial)
    plt.hlines(0, 0, len(displacements_radial)-1, linestyles="dashed")
    plt.ylim([-1.2,1.2])
    plt.savefig(out_file_radial)

  def save_displacements_csv(self, displacementX_file, displacementY_file, out_file_long, out_file_radial, cacheLoc = None):
    displacements_long, displacements_radial = \
      self._load_displacements_through_curve(
        displacementX_file=displacementX_file,
        displacementY_file=displacementY_file,
        cacheLoc=cacheLoc
      )
    
    write_csv_file(displacements_long, out_file_long)
    write_csv_file(displacements_radial, out_file_radial)

  

  
  

  
