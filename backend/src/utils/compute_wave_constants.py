import numpy as np
from scipy.stats import linregress

def _find_indices(array):
  return np.argwhere(array).T[0]

def get_inliers(x, y):
  first_slope, first_intercept, _, _, _ = linregress(x, y)
  first_linear_y = first_slope * x + first_intercept

  errors = np.abs(first_linear_y - y)
  std_dev = np.sqrt(np.sum(errors ** 2) / len(errors))

  not_noisy_indices = _find_indices(errors <= std_dev + 1e-6)
  
  return not_noisy_indices

def _robust_linear_fit(x, y):
  assert x.shape == y.shape
  inliers = get_inliers(x, y)
  final_slope, final_intercept, _, _, _ = linregress(x[inliers], y[inliers])
  return final_slope, final_intercept

def _get_weighted_average(y, weights):
  return np.sum(y * weights) / weights.sum()

def _robust_constant_fit(y, weights):
  weighted_average = _get_weighted_average(y, weights)
  errors = np.abs(y - weighted_average)
  std_dev = np.sqrt(_get_weighted_average(errors**2, weights))
  inliers = _find_indices(errors <= std_dev + 1e-6)
  return _get_weighted_average(y[inliers], weights[inliers])

def compute_difference(phases):
  out = phases[1:] - phases[:-1]
  out %= 2*np.pi
  out[out > np.pi] -= 2*np.pi
  return out

def wave_speed_from_slope(slope, pixel_size, frequency):
  return abs(2 * np.pi * pixel_size * frequency / slope) * 1e-6

def decay_constant_from_slope(slope, pixel_size):
  return -pixel_size / slope

def fit_phases(phases, input_weights):
  weights = input_weights + 1e-6

  first_differences = compute_difference(phases)
  mean_difference = _robust_constant_fit(first_differences, weights[:-1])
  slope = mean_difference

  xs = np.arange(len(phases))
  (sin_intercept, cos_intercept) = weighted_least_squares(
    A = np.array([np.cos(xs * slope), np.sin(xs * slope)]).T,
    b = np.sin(phases),
    weights = weights
  )

  intercept = np.arctan2(sin_intercept, cos_intercept)
  return slope, intercept

def weighted_least_squares(A, b, weights):
  first_x, _, _, _ = np.linalg.lstsq(
    np.dot(np.diag(weights), A),
    weights * b
  )

  estimates = np.dot(A, first_x)
  errors = (b - estimates)
  std_dev = np.sqrt(_get_weighted_average(errors**2, weights))
  inliers = _find_indices(errors <= std_dev + 1e-6)

  final_x, _, _, _ = np.linalg.lstsq(
    np.dot(np.diag(weights), A)[inliers, :],
    (weights * b)[inliers]
  )

  return final_x

def fit_line(y, lower_bound):
  xs = _find_indices(y >= lower_bound)
  if len(xs) < len(y) / 3:
    xs = np.arange(len(y))
  slope, intercept = _robust_linear_fit(xs, y[xs])
  return slope, intercept
