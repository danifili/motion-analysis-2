from src.post_processing.motion_across_curve.frame_to_frame_displacements import DisplacementThroughCurve
from src.motion_analysis.MyImage import MyImage
import sys
import argparse
import os
from tqdm import tqdm


def check_file_exists(filename):
  if not os.path.isfile(filename):
    raise RuntimeError(filename + " does not exists")

def check_dir_exists(root):
  directory = os.path.dirname(root)
  if not os.path.exists(directory):
    raise RuntimeError(directory + " does not exists")

#classes that wrap logic when parsing line in input_file
class MotionAnalysisInput:
  def __init__(self, words):
    assert len(words) in (10, 12)
    self.images = words[:8]
    self.frequency = int(words[8])
    self.out_root = str(words[9])
    self.bounds = tuple(map(int, (words[10:12]))) if len(words) == 12 else None
    #checking existance
    for image in self.images:
      check_file_exists(image)
    check_dir_exists(self.out_root)
  
  def get_script_params(self, args, state):
    flags = []
    #positional arguments
    flags.append(" ".join(self.images))
    flags.append(self.out_root)
    
    #standard flags
    flags.append("-q 0 -x -p --fullImage")
    
    #optional flags
    if not args["fullAmpPhasePlot"]:
      flags.append("-s")
    if not args["fasterAlgorithm"]:
      flags.append("--hornShunck")
    if args["curveDisplacementPlot"] or args["curveDisplacementCSV"]:
      flags.append("--includeDisplacements")

    return " ".join(flags)
    

  def get_curve_post_analysis_input(self):
    words = []
    words.append(self.out_root + "amplitudes_x.csv")
    words.append(self.out_root + "amplitudes_y.csv")
    words.append(self.out_root + "phases_x.csv")
    words.append(self.out_root + "phases_y.csv")
    words.append(self.frequency)
    words.append(self.out_root)
    if self.bounds is not None:
      words.extend(map(str, self.bounds))
    return CurvePostAnalysisInput(words)
  
  def get_curve_displacement_input(self):
    words = []
    for i in range(7):
      words.append(self.out_root + "displacementsX" + str(i) + "to" + str(i+1) + ".csv")
    for i in range(7):
      words.append(self.out_root + "displacementsY" + str(i) + "to" + str(i+1) + ".csv")
    words.append(self.images[0])
    words.append(self.out_root)
    return CurveDisplacementInput(words)

class CurvePostAnalysisInput:
  def __init__(self, words):
    assert len(words) in (6, 8)
    self.ampX = words[0]
    self.ampY = words[1]
    self.phaseX = words[2]
    self.phaseY = words[3]
    self.frequency = int(words[4])
    self.out_root = words[5]
    self.bounds = map(int, (words[6:8])) if len(words) == 8 else None

    #check existance
    # for filename in (self.ampX, self.ampY, self.phaseX, self.phaseY):
    #   check_file_exists(filename)
    check_dir_exists(self.out_root)
  
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
    if self.bounds is not None:
      flags.append("--bounds")
      flags.extend(map(str, self.bounds))
    
    return " ".join(flags)

class CurveDisplacementInput:
  def __init__(self, words):
    assert len(words) == 16
    self.displacementsX = words[:7]
    self.displacementsY = words[7:14]
    self.background_image = MyImage(words[14])
    self.out_root = words[15]

  
#flag handlers
def extract_metadata(args, state):
  with open(args["input_file"], 'r') as f:
    line_number = 0
    state["motion_analysis_input"] = []
    state["curve_post_analysis_input"] = []
    state["curve_displacement_input"] = []
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
          if words[0] == "images":
            motion_inp = MotionAnalysisInput(words[1:])
            state["motion_analysis_input"].append(motion_inp)
            state["curve_post_analysis_input"].append(motion_inp.get_curve_post_analysis_input())
            if args["curveDisplacementPlot"] or args["curveDisplacementCSV"]:
              state["curve_displacement_input"].append(motion_inp.get_curve_displacement_input())
          elif words[0] == "motion":
            post_inp = CurvePostAnalysisInput(words[1:])
            state["curve_post_analysis_input"].append(post_inp)
          elif words[0] == "displacements":
            disp_inp = CurveDisplacementInput(words[1:])
            state["curve_displacement_input"].append(disp_inp)
          else:
            raise RuntimeError("incorrect first argument")
        except:
          raise RuntimeError("Parsing error in line: {}".format(str(line_number)))
      line_number += 1

def analyse_motion(args, state):
  command_prefix = "bash ./scripts/run.sh backend/src/motion_analysis/motion_analysis.py"
  print ("computing overall motion")
  #run processes in parallel
  for inp in tqdm(state["motion_analysis_input"]):
    os.system("{} {}".format(command_prefix, inp.get_script_params(args, state)))
  

def motion_through_curve(args, state):
  command_prefix = "bash ./scripts/run.sh backend/src/post_processing/motion_across_curve/main.py"
  print ("compute motion through curve")
  #run processes in parallel
  for inp in tqdm(state["curve_post_analysis_input"]):
    os.system("{} {}".format(command_prefix, inp.get_script_params(args, state)))

def displacements_through_curve(args, state):
  for inp in tqdm(state["curve_displacement_input"]):
    displacementThroughCurve = DisplacementThroughCurve(
      start_point=state["start_point"],
      end_point=state["end_point"],
      control_point=state["control_point"],
      background_image=inp.background_image,
      pixel_size=state["pixel_size"],
      displacementX_files=inp.displacementsX,
      displacementY_files=inp.displacementsY,
      cacheLoc="backend/cache.pkl"
    )
    for i in range(8):
      if args["curveDisplacementPlot"]:
        displacementThroughCurve.plot_displacements(
          time=i,
          out_file_long=inp.out_root + "displacementsLong" + str(i) + ".eps",
          out_file_radial=inp.out_root + "displacementsRadial" + str(i) + ".eps",
        )
    for i in range(7):
      if args["curveDisplacementCSV"]:
        displacementThroughCurve.save_displacements_csv(
          frame_to_frame_index=i,
          out_file_long=inp.out_root + "displacementsLong" + str(i) + "to" + str(i+1) + ".csv",
          out_file_radial=inp.out_root + "displacementsRadial" + str(i) + "to" + str(i+1) + ".csv",
        )

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
      "The results might have a worse quality but will take less time to compute.",
    "curveDisplacementPlot": "If included, the script will output 14 plots representing the longitudinal and radial " + \
      "displacements alongisde the selected curve",
    "curveDisplacementCSV": "If included, the script will output 14 csv files representing the longitudinal and radial " + \
      "displacements alongside the selected curve. The unit is in pixels."
  }

  flags = [
    {"input_file": dict(action="store", help=HELP["input_file"])},
    {"--curveAmpPhaseCSV": dict(action="store_true", help=HELP["curveAmpPhaseCSV"])},
    {"--curveAmpPhasePlot": dict(action="store_true", help=HELP["curveAmpPhasePlot"])},
    {"--fullAmpPhasePlot": dict(action="store_true", help=HELP["fullAmpPhasePlot"])},
    {"--fasterAlgorithm": dict(action="store_true", help=HELP["fasterAlgorithm"])},
    {"--curveDisplacementPlot": dict(action="store_true", help=HELP["curveDisplacementPlot"])},
    {"--curveDisplacementCSV": dict(action="store_true", help=HELP["curveDisplacementCSV"])}
  ]

  functions = [extract_metadata, analyse_motion, motion_through_curve, displacements_through_curve]

  parser = argparse.ArgumentParser()

  for flag in flags:
      key = list(flag.keys())[0]
      parser.add_argument(key, **flag[key])

  args = vars(parser.parse_args(sys.argv[1:]))
  state = {}
  
  for f in functions:
      f(args, state)