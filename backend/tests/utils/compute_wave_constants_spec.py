import unittest
import numpy as np

from src.utils.compute_wave_constants import *

class TestUnwrappPhase(unittest.TestCase):
  def _compare_arrays(self, a, b):
    assert a.shape == b.shape
    for i in range(len(a)):
      self.assertAlmostEqual(a[i], b[i], 3)
    
  def mod(self, a, b):
    out = a % b
    if out < 0:
      out += b
    return out

  def test_no_unwrapping(self):
    phases = 2 * np.pi * np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    expected = phases
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)
  
  def test_normalizing_input(self):
    phases = 2 * np.pi * np.array([1.2, 2.2, -.8, -1.8])
    expected = 2 * np.pi * np.array([0.2, 0.2, 0.2, 0.2])
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)
  
  def test_basic_upwards_wrapping(self):
    phases = 2 * np.pi * np.array([0.6, 0.8, 0.1, 0.2, 0.4])
    expected = 2 * np.pi * np.array([0.6, 0.8, 1.1, 1.2, 1.4])
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)

  def test_basic_downwards_wrapping(self):
    phases = 2 * np.pi * np.array([0.5, 0.2, 0.7, 0.4, 0.1])
    expected = 2 * np.pi * np.array([0.5, 0.2, -0.3, -0.6, -0.9])
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)
  
  def test_zigzag(self):
    phases = 2 * np.pi * np.array([0.8, 0.2, 0.9])
    expected = 2 * np.pi * np.array([0.8, 1.2, 0.9])
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)
  
  def test_big_upwards_unwrapping(self):
    phases = 2 * np.pi * (0.3 - (np.arange(100) * 0.2) % (1.0))
    expected = 2 * np.pi * (0.3 - (np.arange(100) * 0.2))
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)

  def test_big_downwards_unwrapping(self):
    phases = 2 * np.pi * (0.3 + (np.arange(100) * 0.2) % (1.0))
    expected = 2 * np.pi * (0.3 + (np.arange(100) * 0.2))
    out = unwrap_phases(phases)
    self._compare_arrays(expected, out)


class TestWaveConstants(unittest.TestCase):
  def test_wave_decay(self):
    pixel_size = 0.3
    amplitudes = 1.5 * np.exp(-0.015 * np.arange(100))
    expected = 1 / 0.05
    out = decay_constant(amplitudes, pixel_size)
    self.assertAlmostEqual(expected, out, 3)
  
  def test_wave_speed(self):
    pixel_size = 0.1
    frequency = 20
    micron_in_meters = 1e-6
    phases = 2 * np.pi * (0.3 - (np.arange(100) * 0.2) % (1.0))
    expected = micron_in_meters * (2 * np.pi) * frequency * pixel_size / (0.2 * (2*np.pi))
    out = wave_speed(phases, pixel_size, frequency)
    self.assertAlmostEqual(expected, out, 3)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)