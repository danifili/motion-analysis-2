import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

from numba import njit

from src.utils.compute_wave_constants import decay_constant_from_slope, wave_speed_from_slope
from src.utils.file_reader import read_csv_file, dict_to_json, write_csv_file
from src.post_processing.motion_across_curve.curve_utils import clamp, get_range, get_points_from_GUI, change_basis, LOWER_BOUND_LOG_AMP, compute_phase_weights, compute_amp_phase_fit


@njit
#functions for computing amp and phases throughout the curve
def get_radial_amp_and_phase(long_vec, amp_vec, phase_vec):
  Ax, Bx = amp_vec[0] * np.cos(phase_vec[0]), amp_vec[0] * np.sin(phase_vec[0])
  Ay, By = amp_vec[1] * np.cos(phase_vec[1]), amp_vec[1] * np.sin(phase_vec[1])
  
  new_Ax, new_Ay = change_basis(long_vec, Ax, Ay)
  new_Bx, new_By = change_basis(long_vec, Bx, By)

  new_amp_x = np.sqrt(new_Ax ** 2 + new_Bx ** 2)
  new_amp_y = np.sqrt(new_Ay ** 2 + new_By ** 2)

  new_phase_x = np.arctan2(new_Bx, new_Ax) % (2 * np.pi)
  new_phase_y = np.arctan2(new_By, new_Ay) % (2 * np.pi)

  return np.array([new_amp_x, new_amp_y]), np.array([new_phase_x, new_phase_y])

def get_long_radial_from_XY(start_point, end_point, control_point, ampX, ampY, phaseX, phaseY, cacheLoc=None, translation=1):
  assert ampX.shape == ampY.shape == phaseX.shape == phaseY.shape

  width, height = ampX.shape

  points = get_points_from_GUI(
    start_point=start_point,
    control_point=control_point,
    end_point=end_point,
    width=width,
    height=height,
    cacheLoc=cacheLoc
  )
  
  long_amps = np.zeros((len(points)-1))
  long_phases = np.zeros((len(points)-1))
  radial_amps = np.zeros((len(points)-1))
  radial_phases = np.zeros((len(points)-1))
  for i in range(len(points)-1):
    point = points[i]
    next_point = points[i+1]
    discrete_x, discrete_y = clamp(int(round(point[0])), 0, width-1), clamp(int(round(point[1])), 0, height-1)

    (long_amp, radial_amp), (long_phase, radial_phase) = get_radial_amp_and_phase(
      long_vec = next_point - point,
      amp_vec = np.array([ampX[discrete_x, discrete_y], ampY[discrete_x, discrete_y]]),
      phase_vec = np.array([phaseX[discrete_x, discrete_y], phaseY[discrete_x, discrete_y]])
    )

    long_amps[i] = long_amp
    long_phases[i] = long_phase
    radial_amps[i] = radial_amp
    radial_phases[i] = radial_phase

  return np.array(long_amps), np.array(radial_amps), np.array(long_phases), np.array(radial_phases)

#script handlers
def generate_data(args, state):
  start_point = np.array(list(map(float, args["start_point"])))
  end_point = np.array(list(map(float, args["end_point"])))
  control_point = np.array(list(map(float, args["control_point"])))

  ampX = read_csv_file(args["ampX"])
  ampY = read_csv_file(args["ampY"])
  phaseX = read_csv_file(args["phaseX"])
  phaseY = read_csv_file(args["phaseY"])


  long_amps, radial_amps, long_phases, radial_phases = get_long_radial_from_XY(
    start_point=start_point,
    end_point=end_point,
    control_point=control_point,
    ampX=ampX,
    ampY=ampY,
    phaseX=phaseX,
    phaseY=phaseY,
    cacheLoc = "./backend/cache.pkl" if args["cache"] else None
  )

  state["ampLong"] = long_amps
  state["ampRadial"] = radial_amps
  state["phaseLong"] = long_phases
  state["phaseRadial"] = radial_phases

  state["ampLongSlope"], state["ampLongInt"], state["phaseLongSlope"], state["phaseLongInt"] = \
    compute_amp_phase_fit(
      get_range(state["ampLong"], args["bounds"]),
      get_range(state["phaseLong"], args["bounds"])
    )

  state["ampRadialSlope"], state["ampRadialInt"], state["phaseRadialSlope"], state["phaseRadialInt"] = \
    compute_amp_phase_fit(
      get_range(state["ampRadial"], args["bounds"]),
      get_range(state["phaseRadial"], args["bounds"])
    )

def standard_out(args, state):

  decay_long = decay_constant_from_slope(state["ampLongSlope"], args["pixel_size"])
  speed_long = wave_speed_from_slope(state["phaseLongSlope"], args["pixel_size"], args["frequency"])
  decay_radial = decay_constant_from_slope(state["ampRadialSlope"], args["pixel_size"])
  speed_radial = wave_speed_from_slope(state["phaseRadialSlope"], args["pixel_size"], args["frequency"])

  json_out = {
    "decayLong": decay_long,
    "decayRadial": decay_radial,
    "speedLong": speed_long,
    "speedRadial": speed_radial,
  }

  dict_to_json(json_out, args["root"] + "wave_results.json")
  

def exportCSV_handler(args, state):
  header = str(args["pixel_size"])
  if args["exportCSV"]:
    write_csv_file(state["ampLong"], args["root"] + "amplitudes_long.csv", header=header)
    write_csv_file(state["ampRadial"], args["root"] + "amplitudes_radial.csv", header=header)
    write_csv_file(state["phaseLong"], args["root"] + "phase_long.csv", header=header)
    write_csv_file(state["phaseRadial"], args["root"] + "phase_radial.csv", header=header)

def save_plot_if_included(args, suffix):
    if args["savePlots"]:
      plt.savefig(args["root"] + suffix)

def plot_phases(phases, phase_slope, phase_int, bounds):

  phases_fit = phase_slope * np.arange(len(phases)) + phase_int
  phases_fit %= 2*np.pi
  phases_fit[phases_fit < 0] += 2*np.pi

  plt.xlabel("Pixels")
  plt.ylabel("Phase")

  plt.plot(phases, label="Original")
  plt.plot(get_range(np.arange(len(phases)), bounds), get_range(phases_fit, bounds), label="Best Linear Fit")
  plt.ylim([0,2*np.pi])
  plt.legend()

def plot_amplitudes(amps, amps_slope, amps_int, bounds):
  log_amps_fit = amps_slope * np.arange(len(np.log(amps))) + amps_int

  plt.xlabel("Pixels")
  plt.ylabel("Log Amplitudes in Pixels")

  plt.plot(np.log(amps), label="Original")
  plt.plot(get_range(np.arange(len(amps)), bounds), get_range(log_amps_fit, bounds), label="Best Linear Fit")
  plt.legend()

def plots_handler(args, state):
  
  if args["showPlots"] or args["savePlots"]:
    plt.figure("Long Amp")
    plot_amplitudes(state["ampLong"], state["ampLongSlope"], state["ampLongInt"], args["bounds"])
    save_plot_if_included(args, "amplitudes_long.png")

    plt.figure("Radial Amp")
    plot_amplitudes(state["ampRadial"], state["ampRadialSlope"], state["ampRadialInt"], args["bounds"])
    save_plot_if_included(args, "amplitudes_radial.png")

    plt.figure("Long Phase")
    plot_phases(state["phaseLong"], state["phaseLongSlope"], state["phaseLongInt"], args["bounds"])
    save_plot_if_included(args, "phases_long.png")

    plt.figure("Radial Phase")
    plot_phases(state["phaseRadial"], state["phaseRadialSlope"], state["phaseRadialInt"], args["bounds"])
    save_plot_if_included(args, "phases_radial.png")


    if args["showPlots"]:
      plt.show()


if __name__ == "__main__":
  
  HELP = {
    "root": "the root used to store all the files to be saved",
    "frequency": "in HZ",
    "pixel_size": "in microns",
    "start_point": "start point parameter extracted from gui",
    "end_point": "end point parameter extracted from gui",
    "control_point": "control point extracted from gui",
    "ampX": "csv file containing amplitudes in X direction",
    "ampY": "csv file containing amplitudes in Y direction",
    "phaseX": "csv file containing phase in X direction",
    "phaseY": "csv file containing phase in Y direction",
    "--exportCSV": "if included, it outputs 4 csv files whose prefix are [amplitudes_long.csv] " + \
                   "[amplitudes_radial.csv], [phases_long.csv], [phases_radial.csv]",
    "--showPlots": "show plots for radial and longitudinal ampltudes and phases",
    "--savePlots": "save plots for radial and longitudinal amplitudes and phases with the suffixes" + \
                   "[amplitudes_long.png], [amplitudes_radial.png], [phases_long.png], [phases_long.png]"
  }

  flags = [
    {"root": dict(action="store", help=HELP["root"])},
    {"frequency": dict(action="store",  help=HELP["frequency"], type=float)},
    {"pixel_size": dict(action="store",  help=HELP["pixel_size"], type=float)},
    {"start_point": dict(action="store", nargs=2, help=HELP["start_point"], type=float)},
    {"end_point": dict(action="store", nargs=2, help=HELP["end_point"], type=float)},
    {"control_point": dict(action="store", nargs=2, help=HELP["control_point"], type=float)},
    {"ampX": dict(action="store", help=HELP["ampX"], type=str)},
    {"ampY": dict(action="store", help=HELP["ampY"], type=str)},
    {"phaseX": dict(action="store", help=HELP["phaseX"], type=str)},
    {"phaseY": dict(action="store", help=HELP["phaseY"], type=str)},
    {"--exportCSV": dict(action="store_true", help=HELP["--exportCSV"])},
    {"--savePlots": dict(action="store_true", help=HELP["--savePlots"])},
    {"--showPlots": dict(action="store_true", help=HELP["--showPlots"])},
    {"--bounds": dict(action="store", nargs=2, type=int)},
    {"--cache": dict(action="store_true")}
  ]

  functions = [generate_data, standard_out, exportCSV_handler, plots_handler]

  parser = argparse.ArgumentParser()

  for flag in flags:
      key = list(flag.keys())[0]
      parser.add_argument(key, **flag[key])

  args = vars(parser.parse_args(sys.argv[1:]))
  state = {}
  
  for f in functions:
      f(args, state)

    


  



  







