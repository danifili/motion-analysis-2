import sys
import argparse
import os
import subprocess

from src.post_processing.motion_across_curve.bezier_curve import BezierCurve


#classes that wrap logic when parsing line in input_file
class MotionAnalysisInput:
  def __init__(self, words):
    self.images = words[:8]
    self.frequency = int(words[8])
    self.out_root = str(words[9])
  
  def get_script_params(self, args, state):
    flags = []
    #positional arguments
    flags.append(" ".join(self.images))
    flags.append(self.out_root)
    
    #standard flags
    flags.append("-q 0 -x -p")
    
    #optional flags
    if args["fullAmpPhasePlot"]:
      flags.append("-s")
    if not args["fasterAlgorithm"]:
      flags.append("--hornShunck")
    
    return " ".join(flags)
    

  def get_curve_post_analysis_input(self):
    words = []
    words.append(self.out_root + "amplitudes_x.csv")
    words.append(self.out_root + "amplitudes_y.csv")
    words.append(self.out_root + "phases_x.csv")
    words.append(self.out_root + "phases_y.csv")
    words.append(self.frequency)
    words.append(self.out_root)
    return words

class CurvePostAnalysisInput:
  def __init__(self, words):
    self.ampX = words[0]
    self.ampY = words[1]
    self.phaseX = words[2]
    self.phaseY = words[3]
    self.frequency = int(words[4])
    self.out_root = words[5]
  
  def get_script_params(self, args, state):
    flags = []
    #positional arguments
    flags.append(self.out_root)
    flags.append(str(self.frequency))
    flags.append(str(state["pixel_size"]))
    flags.append(" ".join(list(map(str, state["start_point"]))))
    flags.append(" ".join(list(map(str, state["end_point"]))))
    flags.append(" ".join(list(map(str, state["control_point"]))))
    flags.append(self.ampX)
    flags.append(self.ampY)
    flags.append(self.phaseX)
    flags.append(self.phaseY)

    flags.append("--cache")
    #optional arguments
    if args["curveAmpPhaseCSV"]:
      flags.append("--exportCSV")
    if args["curveAmpPhasePlot"]:
      flags.append("--savePlots")
    
    return " ".join(flags)


#flag handlers
def extract_metadata(args, state):
  with open(args["input_file"], 'r') as f:
    line_number = 0
    state["motion_analysis_input"] = []
    state["curve_post_analysis_input"] = []
    for row in f:
      words = row.split()
      if line_number == 0:
        state["pixel_size"] = float(words[0])
        state["start_point"] = tuple(map(float, words[1:3]))
        state["end_point"] = tuple(map(float, words[3:5]))
        state["control_point"] = tuple(map(float, words[5:7]))
        assert state["pixel_size"] > 0
      else:
        try:
          if len(words) == 10:
            motion_inp = MotionAnalysisInput(words)
            state["motion_analysis_input"].append(motion_inp)
            state["curve_post_analysis_input"].append(motion_inp.get_curve_post_analysis_input())
          elif len(words) == 6:
            post_inp = CurvePostAnalysisInput(words)
            state["curve_post_analysis_input"].append(post_inp)
          else:
            raise RuntimeError("incorrect number of arguments line")
        except:
          raise RuntimeError("Parsing error in line: {}".format(str(line_number)))
      line_number += 1

def analyse_motion(args, state):
  command_prefix = "bash ./scripts/run.sh backend/src/motion_analysis/motion_analysis.py"
  #run processes in parallel
  for inp in state["motion_analysis_input"]:
    os.system(subprocess.Popen("{} {}".format(command_prefix, inp.get_script_params(args, state))))
  

def motion_through_curve(args, state):
  command_prefix = "bash ./scripts/run.sh backend/src/post_processing/motion_across_curve/main.py"

  #run processes in parallel
  for inp in state["curve_post_analysis_input"]:
    os.system("{} {}".format(command_prefix, inp.get_script_params(args, state)))

if __name__ == "__main__":

  #TODO: add flags cache, downsample, motionmag and motionfactor
  HELP = {
    "input_file": "input file containing batches of images to be processed",
    "curveAmpPhaseCSV" : "If specified, the script will output 4 csv files, whose suffixes will be " + \
      "[phase_radial.csv], [phase_long.csv], [amplitudes_radial.csv], [amplitudes_long.csv]",
    "curveAmpPhasePlot":  "If specified, the script will output 4 plots: " + \
      "[phase_radial.png], [phase_long.png], [amplitudes_radial.png], [amplitudes_long.png]",
    "fullAmpPhasePlot": "If specified and images provided, the script will output 4 2-dimensional heat map plots: " + \
      "[amplitudes_x.png], [amplitudes_y.png], [phases_x.png] and [phases_y.png]",
    "fasterAlgorithm": "If included, it will ignore the last step for computing the displacements (The Horn and Schunk algorithm)." + \
      "The results might have a worse quality but will take less time to compute."
  }

  flags = [
    {"input_file": dict(action="store", help=HELP["input_file"])},
    {"--curveAmpPhaseCSV": dict(action="store_true", help=HELP["curveAmpPhaseCSV"])},
    {"--curveAmpPhasePlot": dict(action="store_true", help=HELP["curveAmpPhasePlot"])},
    {"--fullAmpPhasePlot": dict(action="store_true", help=HELP["fullAmpPhasePlot"])},
    {"--fasterAlgorithm": dict(action="store_true", help=HELP["fasterAlgorithm"])}
  ]

  functions = [extract_metadata, analyse_motion, motion_through_curve]

  parser = argparse.ArgumentParser()

  for flag in flags:
      key = list(flag.keys())[0]
      parser.add_argument(key, **flag[key])

  args = vars(parser.parse_args(sys.argv[1:]))
  state = {}
  
  for f in functions:
      f(args, state)