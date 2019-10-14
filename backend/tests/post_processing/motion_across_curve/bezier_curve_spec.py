import unittest
import numpy as np

from src.post_processing.motion_across_curve.bezier_curve import BezierCurve

class TestBezierCurve(unittest.TestCase):
  #Helper functions
  def _compare_points(self, a, b):
    self.assertAlmostEqual(a[0], b[0], 3)
    self.assertAlmostEqual(a[1], b[1], 3)

  #get_parametrized_points test cases
  def test_get_parametrized_point_start_point(self):
    curve = BezierCurve(
      start_point = np.array([10, 30]),
      end_point = np.array([40, 60]),
      control_point = np.array([100, 20])
    )
    self._compare_points(curve.start_point, curve.get_parametrized_point(0))
  
  def test_get_parametrized_point_end_point(self):
    curve = BezierCurve(
      start_point = np.array([30, 30]),
      end_point = np.array([40, 20]),
      control_point = np.array([10, 10])
    )
    self._compare_points(curve.end_point, curve.get_parametrized_point(1))
  
  def test_get_parametrized_point_control_point(self):
    start = np.array([10, 10])
    end = np.array([80, 20])
    control = (start + end) / 2
    
    curve = BezierCurve(
      start_point = start,
      control_point = control,
      end_point = end
    )

    self._compare_points(curve.control_point, curve.get_parametrized_point(0.5))

  def test_get_parametrized_point_random_points(self):
    curve = BezierCurve(
      start_point = np.array([0, 0]),
      control_point = np.array([50, 50]),
      end_point = np.array([100, 100])
    )

    for t in [0.2, 0.4, 0.6, 0.8]:
      self._compare_points(curve.get_parametrized_point(t), np.array([100 * t, 100 * t]))
  
  #split_curve_equally test cases
  def test_split_curve_equally_vertical_line(self):
    for t in [0.2, 0.4, 0.6, 0.8]:
      curve = BezierCurve(
        start_point = np.array([0, 0]),
        control_point = t * np.array([101, 0]),
        end_point = np.array([101, 0])
      )

      ans = np.array([[i, 0] for i in range(101, 5)])
      res = curve.split_curve_equally(5)
      for idx in range(len(ans)):
        self._compare_points(ans[idx], res[idx])
  
  def test_split_curve_equally_horizontal_line(self):
    for t in [0.2, 0.4, 0.6, 0.8]:
      curve = BezierCurve(
        start_point = np.array([0, 0]),
        control_point = t * np.array([0, 101]),
        end_point = np.array([0, 101])
      )

      ans = np.array([[0, i] for i in range(101, 5)])
      res = curve.split_curve_equally(5)
      for idx in range(len(ans)):
        self._compare_points(ans[idx], res[idx])

  #estimate_curve_length test cases
  def test_estimate_curve_length_diagonal(self):
    for t in [0.2, 0.4, 0.6, 0.8]:
      curve = BezierCurve(
        start_point = np.array([0, 0]),
        control_point = t * np.array([100, 100]),
        end_point = np.array([100, 100])
      )

      self.assertAlmostEqual(curve.estimate_curve_length(), 100 * np.sqrt(2), 3)
  
  def test_estimate_curve_length_vertical(self):
    for t in [0.2, 0.4, 0.6, 0.8]:
      curve = BezierCurve(
        start_point = np.array([0, 0]),
        control_point = t * np.array([0, 100]),
        end_point = np.array([0, 100])
      )

      self.assertAlmostEqual(curve.estimate_curve_length(), 100, 3)

  def test_estimate_curve_length_horizontal(self):
    for t in [0.2, 0.4, 0.6, 0.8]:
      curve = BezierCurve(
        start_point = np.array([0, 0]),
        control_point = t * np.array([100, 0]),
        end_point = np.array([100, 0])
      )

      self.assertAlmostEqual(curve.estimate_curve_length(), 100, 3)
  
  def test_estimate_curve_length_not_mid_points(self):
    for t1 in [0.2, 0.4]:
      for t2 in [0.6, 0.8]:
        curve = BezierCurve(
          start_point = np.array([0, 0]),
          control_point = np.array([50, 50]),
          end_point = np.array([100, 100])
        )
        self.assertAlmostEqual(curve.estimate_curve_length(start_t=t1, end_t=t2), (t2-t1) * 100 * np.sqrt(2), 3)
  
  #find_parameter_t test cases
  def test_find_parameter_t_diagonal(self):
    curve = BezierCurve(
      start_point = np.array([0, 0]),
      control_point = np.array([50, 50]),
      end_point = np.array([100, 100])
    )

    for t in [0.2, 0.4, 0.6, 0.8]:
      self.assertAlmostEqual(curve.find_parameter_t(t * np.sqrt(2) * 100), t, 3)
  
  def test_find_parameter_t_integration(self):
    curve = BezierCurve(
      start_point = np.array([10, 30]),
      end_point = np.array([40, 60]),
      control_point = np.array([100, 20])
    )

    for t in [0.2, 0.4, 0.6, 0.8]:
      start_t = 0.1
      end_t = 0.1 + t
      self.assertAlmostEqual(
        curve.find_parameter_t(curve.estimate_curve_length(start_t = start_t, end_t = end_t), start_t = start_t),
        end_t,
        places = 3
      )

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
