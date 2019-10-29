import numpy as np

class BezierCurve:
  """
  2D-Quadratic bezier curve determined by three points.
  This curve is parametrized with the following equation.
  P(t) = (1-t)**2 * P0 + [2 * t * (1-t)] * P1 + t ** 2 * P2
  where 0 <= t <= 1 and P0, P1 and P2 are the start, control
  and endpoint respectively
  """
  def __init__(self, start_point, end_point, control_point):
    self.start_point = np.array(start_point)
    self.end_point = np.array(end_point)
    self.control_point = np.array(control_point)
  
  def split_curve_equally(self, curve_segment_length=1):
    """
    Splits curve in chucks of (approximate) equal lengths. The length of the 
    chuncks is determined by the [curve_segment_length] param
    """
    curve_length = self.estimate_curve_length()
    num_steps = int(curve_length / curve_segment_length)
    out = np.array([self.get_parametrized_point(self.find_parameter_t(curve_length * step/num_steps)) for step in range(num_steps)])
    return out
    
  def estimate_curve_length(self, start_t=0, end_t=1, num_steps=100):
    """
    Estimates sub-curve's length from P([start_t]) to P([end_t]) by taking [num_steps] points in the sub curve
    """
    points = [self.get_parametrized_point(start_t + (end_t - start_t)*step/num_steps) for step in range(num_steps+1)]
    distance = 0
    for idx in range(len(points)-1):
      distance += np.linalg.norm(points[idx] - points[idx+1])
    return distance
    
  def get_parametrized_point(self, t):
    return (1-t)**2 * self.start_point + 2*(1-t)*t*self.control_point + t**2* self.end_point
  
  def find_parameter_t(self, target_distance, start_t=0, end_t=1, eps = 1e-6):
    """
    Finds parameter t that corresponds to the point P(t) in sub-curve between P([start_t]) and P([end_t])
    that is target_distance away from P([start_t]). It does so by using binary search. The parameter
    [eps] determines how quickly the search ends. Lower [eps] gives better answers but takes longer.
    """
    mid_t = (start_t + end_t) / 2
    distance_to_mid_t = self.estimate_curve_length(start_t=start_t, end_t=mid_t)
    offset = target_distance - distance_to_mid_t
    if np.abs(offset) < eps:
      return mid_t
    elif offset < 0:
      return self.find_parameter_t(target_distance, start_t=start_t, end_t=mid_t, eps=eps)
    else:
      return self.find_parameter_t(target_distance-distance_to_mid_t, start_t=mid_t, end_t=end_t, eps=eps)