import numpy as np

def decay_constant(amplitudes, pixel_size):
  """
  Computes decay constant

  Params:
    amplitudes: in pixels
    pixel_size: in microns
  
  Returns: 
    decay constant in 1/microns
  """
  log_amps = np.log(amplitudes)
  slope, offset = np.polyfit(np.arange(len(amplitudes)), log_amps, deg=1)
  return -pixel_size / slope

def unwrap_phases(raw_phases):
  assert len(raw_phases) > 0, "input phases must not be empty"
  phases = np.array(raw_phases)
  phases %= 2 * np.pi
  phases[phases < 0] += 2 * np.pi
  out = []
  for i in range(len(phases)):
    if i == 0:
      out.append(phases[0])
    else:
      previous_phase = out[i-1]
      previous_phase_mod_2pi = phases[i-1]
      previous_bucket_2pi = (previous_phase - previous_phase_mod_2pi) / (2 * np.pi)
      new_phase_candidates = [(previous_bucket_2pi + delta) * (2*np.pi) + phases[i] for delta in (-1,0,1)]
      new_phase = min(new_phase_candidates, key=lambda x: abs(x-previous_phase))
      out.append(new_phase)
  return np.array(out)

def wave_speed(phases, pixel_size, frequency):
  """
  Computes wave speed in m/s

  Params:
    phases: values between 0 and 2pi
    pixel_size: in microns
    frequency: in Hz
  
  Returns:
    positive float in m/s
  """
  unwrapped_phases = unwrap_phases(phases)
  slope, offset = np.polyfit(np.arange(len(unwrapped_phases)), unwrapped_phases,  deg=1)
  return abs(2*np.pi * pixel_size * frequency / slope) * 1e-6
  

def fit_line(y):
  slope, offset = np.polyfit(np.arange(len(y)), y, 1)
  if slope == 0:
    slope = 1e-6
  return slope * np.arange(len(y)) + offset

