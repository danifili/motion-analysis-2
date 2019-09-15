import numpy as np

class BezierCurve:
  def __init__(self, start_point, end_point, control_point):
    self.start_point = np.array(start_point)
    self.end_point = np.array(end_point)
    self.control_point = np.array(control_point)
  
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
    return (1-t)**2 * self.start_point + 2*(1-t)*t*self.control_point + t**2* self.end_point
  
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