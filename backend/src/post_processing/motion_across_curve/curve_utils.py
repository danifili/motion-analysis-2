from src.post_processing.motion_across_curve.bezier_curve import BezierCurve
from src.utils.compute_wave_constants import fit_line, fit_phases
from scipy.optimize import least_squares

import numpy as np
from numba import njit

LOWER_BOUND_LOG_AMP = -2

def get_range(array, bounds):
  if bounds is None:
    return array
  min_x, max_x = bounds
  assert 0 <= min_x < len(array)
  assert 0 <= min_x < max_x <= len(array)
  return array[min_x: max_x]

def clamp(x, lower, upper):
  if x < lower:
    return lower
  if x > upper:
    return upper
  return x

def get_points_from_GUI(start_point, end_point, control_point, width, height, cacheLoc):
  real_start_point = start_point * np.array([width/100, height/100])
  real_control_point = control_point * np.array([width/100, height/100])
  real_end_point = end_point * np.array([width/100, height/100])
  curve_segment_length = 1

  curve = BezierCurve(
    start_point = real_start_point,
    control_point = real_control_point,
    end_point = real_end_point
  )

  if cacheLoc is not None:
    return curve.split_curve_equally_with_cache(curve_segment_length, cacheLoc)
  else:
    return curve.split_curve_equally(curve_segment_length)

@njit
def change_basis(reference_long_vec, vecX, vecY):
  theta = np.arctan2(reference_long_vec[1], reference_long_vec[0])
  long_vec = vecX * np.cos(theta) + vecY * np.sin(theta)
  radial_vec = -vecX * np.sin(theta) + vecY * np.cos(theta)
  return long_vec, radial_vec

def compute_phase_weights(amps):
  weights = np.log(amps) - LOWER_BOUND_LOG_AMP
  weights[weights < 0] = 0
  return weights + 1e-6

def compute_amp_phase_fit(amps, phases):
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
