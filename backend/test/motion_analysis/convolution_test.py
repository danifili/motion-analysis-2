import unittest
import numpy as np

from src.motion_analysis.MyVideoHelper2 import optimized_conv

class TestConvolution(unittest.TestCase):
  def test_01(self):
    res = optimized_conv(
      A = np.array([[
        [0.0] * 7,
        [1.0] * 7,
        [2.0] * 7,
        [3.0] * 7,
        [4.0] * 7,
        [5.0] * 7,
        [6.0] * 7,
      ] for i in range(8)]),
      kernel = np.array(
        [
          [0, 1, 0],
          [0, 0, 0],
          [0, 0, 0]
        ]
      ),
      last_kernel = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    )
    ans = np.array([
        [3.0] * 3,
        [4.0] * 3,
        [5.0] * 3,
    ])

    self.assertEqual(res.shape, ans.shape)
    for x in range(res.shape[0]):
      for y in range(res.shape[1]):
        self.assertEqual(res[x, y], ans[x, y])

if __name__ == '__main__':
  res = unittest.main(verbosity=3, exit=False)


