from src.utils.file_reader import read_csv_file, write_csv_file
from src.post_processing.motion_across_curve.bezier_curve import BezierCurve
from src.post_processing.motion_across_curve.curve_utils import *

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
class DisplacementThroughCurve:
  def __init__(self, start_point, end_point, control_point, background_image, pixel_size):
    self.start_point = start_point
    self.end_point = end_point
    self.control_point = control_point
    self.background_image = background_image
    self.pixel_size = pixel_size
  
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

    width, height = read_csv_file(displacementX_file).shape

    plt.figure()

    plt.xlabel("Micrometers")
    plt.ylabel("Displacement in Micrometers")
    self._graph_plot(displacements_long, width)
    plt.savefig(out_file_long, format="eps")

    plt.figure()
    plt.xlabel("Micrometers")
    plt.ylabel("Displacement in Micrometers")
    self._graph_plot(displacements_radial, width)
    plt.savefig(out_file_radial, format="eps")

  def _to_plot_bounds(self, xs, image_width, curve_length):
    return self.start_point[0] / 100 * image_width + (xs / curve_length * (self.end_point[0] - self.start_point[0]) / 100 * image_width)

  def _to_micrometers(self, xs):
    return xs * self.pixel_size

  def _graph_plot(self, y, image_width):
    sample_size = int(0.03 * self._to_plot_bounds(len(y), image_width, len(y))) + 1
    sample_xs = np.arange(0, len(y)-1, sample_size)
    sample_ys = y[::sample_size]
    
    #smooth fit
    xs = np.arange(len(y))
    params = np.polyfit(xs, y, deg=5)
    plt.plot(self._to_micrometers(self._to_plot_bounds(xs, image_width, len(y))), self._to_micrometers(sum(params[5-i] * xs ** i for i in range(6))))

    #scatter plot
    plt.scatter(self._to_micrometers(self._to_plot_bounds(sample_xs, image_width, len(y))), self._to_micrometers(sample_ys), facecolors='none', edgecolors="b")
    plt.hlines(0, 0, image_width-1, linestyles="dashed")

    #draw background image
    ext = np.array([0, image_width, -1.2, 1.2])
    plt.imshow(self.background_image.get_pixels().T, cmap='gray',vmin=0,vmax=255, extent=self._to_micrometers(ext))
    aspect=self.background_image.height/float(self.background_image.width)*((ext[1]-ext[0])/(ext[3]-ext[2]))
    plt.gca().set_aspect(aspect)

  def save_displacements_csv(self, displacementX_file, displacementY_file, out_file_long, out_file_radial, cacheLoc = None):
    displacements_long, displacements_radial = \
      self._load_displacements_through_curve(
        displacementX_file=displacementX_file,
        displacementY_file=displacementY_file,
        cacheLoc=cacheLoc
      )
    
    write_csv_file(displacements_long, out_file_long)
    write_csv_file(displacements_radial, out_file_radial)

  

  
  

  
