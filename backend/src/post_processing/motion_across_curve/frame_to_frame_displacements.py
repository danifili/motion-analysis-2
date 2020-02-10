from src.utils.file_reader import read_csv_file, write_csv_file
from src.post_processing.motion_across_curve.curve_utils import *
from src.motion_analysis.MyVideoFinal import sinusoidal_fit
from src.utils.compute_wave_constants import fit_line, fit_phases
from scipy.optimize import least_squares
import numpy as np
import matplotlib.pyplot as plt
class DisplacementThroughCurve:
  def __init__(
      self, 
      start_point, 
      end_point, 
      control_point, 
      background_image, 
      pixel_size,
      displacementX_files,
      displacementY_files,
      cacheLoc = None
    ):
    self.start_point = start_point
    self.end_point = end_point
    self.control_point = control_point
    self.background_image = background_image
    self.pixel_size = pixel_size

    self.displacements_long = []
    self.displacements_radial = []
    
    num_consecutive_frames = 7
    for i in range(num_consecutive_frames):
      long_dis, radial_dis = self._load_displacements_through_curve(
        displacementX_file = displacementX_files[i],
        displacementY_file = displacementY_files[i],
        cacheLoc = cacheLoc
      )
      self.displacements_long.append(long_dis)
      self.displacements_radial.append(radial_dis)
    
    self.displacements_long = np.array(self.displacements_long)
    self.displacements_radial = np.array(self.displacements_radial)

    self.amp_long, self.amp_radial, self.phase_long, self.phase_radial = self._compute_amp_and_phase()

    self.amp_long_slope, self.amp_long_intercept, self.phase_long_slope, self.phase_long_intercept = \
      self._compute_amp_phase_fit(self.amp_long, self.phase_long)

    self.amp_radial_slope, self.amp_radial_intercept, self.phase_radial_slope, self.phase_radial_intercept = \
      self._compute_amp_phase_fit(self.amp_radial, self.phase_radial)

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
  
  def plot_displacements(self, time, out_file_long, out_file_radial):
    cumulative_displacements = self._cumulative_displacements()

    displacements_long = cumulative_displacements[time, :, 0]
    displacements_radial = cumulative_displacements[time, :, 1]

    plt.figure()
    plt.xlabel("Micrometers")
    plt.ylabel("Displacement in Micrometers")
    self._graph_plot(
      displacements_long,
      self.amp_long_slope,
      self.phase_long_slope,
      self.amp_long_intercept,
      self.phase_long_intercept,
      self.amp_long,
      self.phase_long,
      time
    )
    plt.savefig(out_file_long, format="eps")

    plt.figure()
    plt.xlabel("Micrometers")
    plt.ylabel("Displacement in Micrometers")
    self._graph_plot(
      displacements_radial, 
      self.amp_radial_slope,
      self.phase_radial_slope,
      self.amp_radial_intercept,
      self.phase_radial_intercept,
      self.amp_radial,
      self.phase_radial,
      time
    )
    plt.savefig(out_file_radial, format="eps")
  
  def _cumulative_displacements(self):
    frames = 8
    curve_length = self.displacements_long.shape[1]

    cumulative_displacements = np.zeros((frames, curve_length, 2))
    displacements = np.concatenate(
      [
        self.displacements_long[:,:,np.newaxis], 
        self.displacements_radial[:,:,np.newaxis]
      ], 
      axis=2
    )
    for i in range(1,frames):
      cumulative_displacements[i] =  cumulative_displacements[i-1] + displacements[i-1]
    
    return np.array(cumulative_displacements)

  def _compute_amp_and_phase(self):
    data_fit = sinusoidal_fit(self._cumulative_displacements()[:,:,np.newaxis,:])
    return data_fit[:, 0, 0], data_fit[:, 0, 1], data_fit[:, 0, 2], data_fit[:, 0, 3]

  def _to_plot_bounds(self, xs, curve_length):
    return self.start_point[0] / 100 * self.background_image.width + (xs / curve_length * (self.end_point[0] - self.start_point[0]) / 100 * self.background_image.width)

  def _to_micrometers(self, xs):
    return xs * self.pixel_size

  def _graph_plot(self, y, amp_slope, phase_slope, amp_intercept, phase_intercept, amps, phases, time):
    sample_size = int(0.03 * self._to_plot_bounds(len(y), len(y))) + 1
    sample_xs = np.arange(0, len(y)-1, sample_size)
    sample_ys = y[::sample_size]

    w = 2 * np.pi / 8

    fitted_log_amps = amp_slope * np.arange(len(y)) + amp_intercept
    fitted_phases = phase_slope * np.arange(len(y)) + phase_intercept

    #smooth fit
    plt.plot(
      self._to_micrometers(self._to_plot_bounds(np.arange(len(y)), len(y))),
      self._to_micrometers(
        np.exp(fitted_log_amps) * np.sin(w*time + fitted_phases)
      )
    )

    #cover
    plt.plot(
      self._to_micrometers(self._to_plot_bounds(np.arange(len(y)), len(y))),
      self._to_micrometers(np.exp(fitted_log_amps)),
      color = "blue",
      linestyle = "dashed",
      linewidth=1
    )

    plt.plot(
      self._to_micrometers(self._to_plot_bounds(np.arange(len(y)), len(y))),
      self._to_micrometers(-np.exp(amp_intercept + amp_slope * np.arange(len(y)))),
      color = "blue",
      linestyle = "dashed",
      linewidth = 1
    )

    #scatter plot
    plt.scatter(self._to_micrometers(self._to_plot_bounds(sample_xs, len(y))), self._to_micrometers(sample_ys + amps[sample_xs] * np.sin(phases[sample_xs])), facecolors='none', edgecolors="b")
    plt.hlines(0, 0, self.background_image.width-1, linestyles="dashed")

    #draw background image
    ext = np.array([0, self.background_image.width, -4/3, 4/3])
    plt.imshow(self.background_image.get_pixels().T, cmap='gray',vmin=0,vmax=255, extent=self._to_micrometers(ext))
    aspect=self.background_image.height/float(self.background_image.width)*((ext[1]-ext[0])/(ext[3]-ext[2]))
    plt.gca().set_aspect(aspect)

  def save_displacements_csv(self, frame_to_frame_index, out_file_long, out_file_radial):
    displacements_long = self.displacements_long[frame_to_frame_index]
    displacements_radial = self.displacements_radial[frame_to_frame_index]
    write_csv_file(displacements_long, out_file_long)
    write_csv_file(displacements_radial, out_file_radial)

  @staticmethod
  def _compute_amp_phase_fit(amps, phases):
    weights = compute_phase_weights(amps)
    def damped_wave(y):
      amp_slope, amp_intercept, phase_slope, phase_intercept = y
      x = np.arange(len(amps))
      return np.concatenate(
        [
          np.sqrt(weights) * (amps * np.sin(phases) - np.exp(amp_intercept + amp_slope * x) * np.sin(phase_slope * x + phase_intercept)),
          np.sqrt(weights) * (amps * np.cos(phases) - np.exp(amp_intercept + amp_slope * x) * np.cos(phase_slope * x + phase_intercept))
        ]
      )

    amp_slope_guess, amp_intercept_guess = fit_line(np.log(amps), LOWER_BOUND_LOG_AMP)
    phase_slope_guess, phase_intercept_guess = fit_phases(phases, weights)

    params = np.array([amp_slope_guess, amp_intercept_guess, phase_slope_guess, phase_intercept_guess])

    return least_squares(damped_wave, params).x
  
  

  
