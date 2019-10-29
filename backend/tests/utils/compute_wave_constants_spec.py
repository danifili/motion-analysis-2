import unittest
import numpy as np

from src.utils.compute_wave_constants import *

class TestWaveConstants(unittest.TestCase):
  def test_wave_decay(self):
    pixel_size = 0.3
    amplitudes = 1.5 * np.exp(-0.015 * np.arange(100))
    expected = 1 / 0.05
    out, _ = decay_constant(amplitudes, pixel_size)
    self.assertAlmostEqual(expected, out, 3)
  
  def test_wave_speed(self):
    pixel_size = 0.1
    frequency = 20
    micron_in_meters = 1e-6
    phases = 2 * np.pi * (0.3 - (np.arange(100) * 0.2) % (1.0))
    expected = micron_in_meters * (2 * np.pi) * frequency * pixel_size / (0.2 * (2*np.pi))
    out, _ = wave_speed(phases, pixel_size, frequency, np.ones(len(phases)))
    self.assertAlmostEqual(expected, out, 3)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)