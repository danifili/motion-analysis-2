from src.post_processing.motion_across_curve.bezier_curve import BezierCurve

import numpy as np
from numba import njit

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

