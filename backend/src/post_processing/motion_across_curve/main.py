import numpy as np
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse

from src.post_processing.motion_across_curve.bezier_curve import BezierCurve
from src.utils.compute_wave_constants import decay_constant, wave_speed, unwrap_phases, fit_line
from src.utils.file_reader import read_csv_file, dict_to_json, write_csv_file

#functions for computing amp and phases throughout the curve
def get_radial_amp_and_phase(long_vec, amp_vec, phase_vec):
  theta = np.arctan2(long_vec[1], long_vec[0])
  Ax, Bx = amp_vec[0] * np.cos(phase_vec[0]), amp_vec[0] * np.sin(phase_vec[0])
  Ay, By = amp_vec[1] * np.cos(phase_vec[1]), amp_vec[1] * np.sin(phase_vec[1])
  
  new_Ax = Ax * np.cos(theta) + Ay * np.sin(theta)
  new_Bx = Bx * np.cos(theta) + By * np.sin(theta)
  new_Ay = -Ax * np.sin(theta) + Ay * np.cos(theta)
  new_By = -Bx * np.sin(theta) + By * np.cos(theta)

  new_amp_x = np.sqrt(new_Ax ** 2 + new_Bx ** 2)
  new_amp_y = np.sqrt(new_Ay ** 2 + new_By ** 2)

  new_phase_x = np.arctan2(new_Bx, new_Ax) % (2 * np.pi)
  new_phase_y = np.arctan2(new_By, new_Ay) % (2 * np.pi)

  return np.array([new_amp_x, new_amp_y]), np.array([new_phase_x, new_phase_y])

def get_long_radial_from_XY(start_point, end_point, control_point, ampX, ampY, phaseX, phaseY):
  assert ampX.shape == ampY.shape == phaseX.shape == phaseY.shape

  width, height = ampX.shape

  bezier_curve = BezierCurve(
    start_point = start_point * np.array([width/100, height/100]),
    control_point = control_point * np.array([width/100, height/100]),
    end_point = end_point * np.array([width/100, height/100])
  )

  points = bezier_curve.split_curve_equally(curve_segment_length=1)

  print("divided bezier curve")

  long_amps = []
  long_phases = []
  radial_amps = []
  radial_phases = []

  for i in tqdm(range(len(points)-1)):
    point = points[i]
    next_point = points[i+1]
    discrete_x, discrete_y = int(round(point[0])), int(round(point[1]))
    
    (long_amp, radial_amp), (long_phase, radial_phase) = get_radial_amp_and_phase(
      long_vec = next_point - point,
      amp_vec = np.array([ampX[discrete_x, discrete_y], ampY[discrete_x, discrete_y]]),
      phase_vec = np.array([phaseX[discrete_x, discrete_y], phaseY[discrete_x, discrete_y]])
    )

    long_amps.append(long_amp)
    long_phases.append(long_phase)
    radial_amps.append(radial_amp)
    radial_phases.append(radial_phase)
  
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

  print ("read csv files")

  long_amps, radial_amps, long_phases, radial_phases = get_long_radial_from_XY(
    start_point=start_point,
    end_point=end_point,
    control_point=control_point,
    ampX=ampX,
    ampY=ampY,
    phaseX=phaseX,
    phaseY=phaseY
  )

  state["ampLong"] = long_amps
  state["ampRadial"] = radial_amps
  state["phaseLong"] = long_phases
  state["phaseRadial"] = radial_phases

def standard_out(args, state):
  decay_long = decay_constant(state["ampLong"], args["pixel_size"])
  decay_radial = decay_constant(state["ampRadial"], args["pixel_size"])
  speed_long = wave_speed(state["phaseLong"], args["pixel_size"], args["frequency"])
  speed_radial = wave_speed(state["phaseRadial"], args["pixel_size"], args["frequency"])

  json_out = {
    "decayLong": decay_long,
    "decayRadial": decay_radial,
    "speedLong": speed_long,
    "speedRadial": speed_radial
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

def plots_handler(args, state):
  
  if args["showPlots"] or args["savePlots"]:
    plt.figure("Long Amp")
    plt.plot(np.log(state["ampLong"]))
    plt.plot(fit_line(np.log(state["ampLong"])))
    save_plot_if_included(args, "amplitudes_long.png")

    plt.figure("Long Phase")
    plt.plot(state["phaseLong"])
    plt.plot(fit_line(unwrap_phases(state["phaseLong"])))
    save_plot_if_included(args, "phases_long.png")

    plt.figure("Radial Amp")
    plt.plot(np.log(state["ampRadial"]))
    plt.plot(fit_line(np.log(state["ampRadial"])))
    save_plot_if_included(args, "amplitudes_radial.png")

    plt.figure("Radial Phase")
    plt.plot(state["phaseRadial"])
    plt.plot(fit_line(unwrap_phases(state["phaseRadial"])))
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
    {"--showPlots": dict(action="store_true", help=HELP["--showPlots"])}
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

    


  



  







