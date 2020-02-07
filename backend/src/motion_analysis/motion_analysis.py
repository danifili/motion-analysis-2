#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 12:04:33 2017

@author: danifili
"""

import sys
from src.motion_analysis.MyVideoFinal import MyVideo, sinusoidal_fit
from src.motion_analysis.MyImage import MyImage
from src.motion_analysis.Plot import Plot
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

METAVAR_IMAGE = tuple("frame" + str(t) for t in range(8))
#helper function
def get_displacements_frame_to_frame(cumulative_displacements):
    duration, width, height, directions = cumulative_displacements.shape
    displacements = np.zeros((duration-1, width, height, directions))
    for i in range(duration-1):
        displacements[i] = cumulative_displacements[i+1] - cumulative_displacements[i]
    return displacements

#action functions
def generate_data(args):
    images = [MyImage(image, compress_image=not args["fullImage"]) for image in args["image"]]
    video = MyVideo(images)

    min_corner = tuple(args["c"][:2])
    x_max, y_max = tuple(args["c"][2:])
    if x_max <= 0:
        x_max += video.width-1
    if y_max <= 0:
        y_max += video.height-1
    max_corner = (x_max, y_max)

    quality_level = args["q"]

    max_iterations = 100 if args["hornShunck"] else 0

    cumulative_displacements = video.get_cumulative_displacements(min_corner, max_corner, quality_level=quality_level, max_iterations=max_iterations)
    data = sinusoidal_fit(cumulative_displacements)
    
    args["data"] = data
    args["cumulative_displacements"] = cumulative_displacements
    args["video"] = video
    args["min_corner"] = min_corner
    args["max_corner"] = max_corner

def save_data(args):
    """
    Saves the amplitudes and phases in the x and y directions of a certain ROI.        
    """
    data = args["data"]
    root = args["root"]
    save_csv = args["x"]

    if not save_csv:
        return

    wildcards = ["amplitudes_x", "amplitudes_y", "phases_x", "phases_y"]

    for i in range(len(wildcards)):
        np.savetxt(root + wildcards[i] + ".csv", data[:, :, i], delimiter=",")

    if args["includeErrors"]:
        np.savetxt(root + "errors_x" + ".csv", data[:, :, 4], delimiter=",")
        np.savetxt(root + "errors_y" + ".csv", data[:, :, 5], delimiter=",")
    
    if args["includeDisplacements"]:
        displacements_frame_to_frame = get_displacements_frame_to_frame(args["cumulative_displacements"])
        assert displacements_frame_to_frame.shape[0] == 7
        assert displacements_frame_to_frame.shape[1:3] == data.shape[:2]
        for i in range(7):
            np.savetxt(root + "displacementsX" + str(i) + "to" + str(i+1) + ".csv", displacements_frame_to_frame[i, :, :, 0], delimiter=",")
            np.savetxt(root + "displacementsY" + str(i) + "to" + str(i+1) + ".csv", displacements_frame_to_frame[i, :, :, 1], delimiter=",")

def get_pixels(video):
    return video[0:video.width, 0:video.height, 0]

def plot_data(args):
    root = args["root"]
    video = args["video"]
    data = args["data"]
    min_corner = args["min_corner"]
    max_corner = args["max_corner"]
    thr_x, thr_y = args["t"]

    arrows_plot = args["a"]
    save_plots = args["s"]
    show_plots = args["p"]
    wave = args["w"] is not None
    arrows_8 = args["a8"]

    cumulative_displacements = args["cumulative_displacements"]

    amplitudes_x = np.array(data[:, :, 0])
    amplitudes_y = np.array(data[:, :, 1])
    phases_x = np.array(data[:, :, 2])
    phases_y = np.array(data[:, :, 3])

    phases_x[amplitudes_x < thr_x] = float('nan')
    amplitudes_x[amplitudes_x < thr_x] = float('nan')

    phases_y[amplitudes_y < thr_y] = float('nan')
    amplitudes_y[amplitudes_y < thr_y] = float('nan')

    if arrows_plot:
        amplitudes = plt.figure()
        k, scale = args["a"]
        Plot.plot(get_pixels(video), data[:,:,:2], min_corner, max_corner, k=int(k), scale=1/scale, color='red')
    
    else:
        amplitudes_x_fig = plt.figure('Amplitudes X')
        plt.title('Amplitudes X')
        Plot.scalar_heat_map(get_pixels(video), min_corner, max_corner, amplitudes_x, alpha=0.3)

        amplitudes_y_fig = plt.figure('Amplitudes Y')
        plt.title('Amplitudes Y')
        Plot.scalar_heat_map(get_pixels(video), min_corner, max_corner, amplitudes_y, alpha=0.3)

        phases_x_fig = plt.figure('Phases X')
        plt.title('Phases X')
        Plot.phase_heat_map(get_pixels(video), min_corner, max_corner, phases_x, alpha=0.3)

        phases_y_fig = plt.figure('Phases Y')
        plt.title('Phases Y')
        Plot.phase_heat_map(get_pixels(video), min_corner, max_corner, phases_y, alpha=0.3)

        error_x_fig = plt.figure('Error X')
        plt.title('Error X')
        Plot.scalar_heat_map(get_pixels(video), min_corner, max_corner, args["data"][:, :, 4], alpha=0.3)
        
        error_y_fig = plt.figure('Error Y')
        plt.title('Error Y')
        Plot.scalar_heat_map(get_pixels(video), min_corner, max_corner, args["data"][:, :, 5], alpha=0.3)

    if wave:
        wave_data = args["wave_data"]
        frequency = args["frequency"]
        pixel_size = args["pixel_size"]
        wave_width, wave_height = wave_data.shape[:2]

        x = np.arange(wave_width)
        y = np.arange(wave_height)
        Ax = np.log(np.array([np.average(wave_data[:, j, 0]) for j in range(wave_height)]))
        Ay = np.log(np.array([np.average(wave_data[i, :, 1]) for i in range(wave_width)]))
        Px = np.unwrap(wave_data[wave_width//2, :, 2])
        Py = np.unwrap(wave_data[:, wave_height//2, 3])

        for (xi, yi, name, dimension) in [(x, Ay, "decay constant x", "1/nm"), 
                                          (y, Ax, "decay constant y", "1/nm"), 
                                          (x, Py, "wave speed x", "m/s"),
                                          (y, Px, "wave speed y", "m/s")]:
            plt.figure()
            plt.plot(xi , yi)
            a , b = np.polyfit(xi, yi, 1)
            plt.plot(xi, a*xi+b*np.ones((len(xi))))
            if dimension == "m/s":
                a = 2*np.pi*frequency / (1e9 * a / pixel_size)
            else:
                a = a / pixel_size

            print (name + ": " + str(abs(a)) + " " + dimension)

    if save_plots:
        if arrows_plot:
            amplitudes.savefig(root + "amplitudes_vector_field.png")
        else:
            amplitudes_x_fig.savefig(root + "amplitudes_x" + ".png")
            amplitudes_y_fig.savefig(root + "amplitudes_y" + ".png")
            phases_x_fig.savefig(root + "phases_x" + ".png")
            phases_y_fig.savefig(root + "phases_y" + ".png")
            error_x_fig.savefig(root + "error_x" + ".png")
            error_y_fig.savefig(root + "error_y" + ".png")

        if arrows_8:
            k8, scale8 = arrows_8
            displacements = get_displacements_frame_to_frame(cumulative_displacements)
            for time in range(len(displacements)):
                optical_flow = plt.figure()
                Plot.plot(video, displacements[time], min_corner, max_corner, time, k=int(k8), scale=1/scale8, color="red")
                optical_flow.savefig(root + "displacements_vector_field_" + str(time) + ".png")
                plt.close(optical_flow)

    if show_plots:
        plt.show()


def motion_mag(args):
    factor = args["motionmag"]
    if factor is None:
        return

    cumulative_displacements = args["cumulative_displacements"]
    video = args["video"]
    root = args["root"]
    x_min, y_min = args["min_corner"]
    x_max, y_max = args["max_corner"]

    u = lambda x, y, t: factor * cumulative_displacements[t, x, y, 0]
    v = lambda x, y, t: factor * cumulative_displacements[t, x, y, 1]

    original_image = MyImage.image_from_matrix(video[x_min:x_max+1, y_min:y_max+1, 0], root + "motion_mag_factor" + str(factor) + "_" + str(0) + ".bmp")

    for t in range(1, video.duration):
        u_t = lambda x, y: u(x, y, t)
        v_t = lambda x, y: v(x, y, t)
        original_image.shift_image(u_t, v_t, root + "motion_mag_factor" + str(factor) + "_" + str(t) + ".bmp")

def motion_stop(args):
    if not args["motionstop"]:
        return

    cumulative_displacements = args["cumulative_displacements"]
    video = args["video"]
    root = args["root"]
    x_min, y_min = args["min_corner"]
    x_max, y_max = args["max_corner"]

    u = lambda x, y, t: -cumulative_displacements[t, x, y, 0]
    v = lambda x, y, t: -cumulative_displacements[t, x, y, 1]
    
    for t in range(video.duration):
        u_t = lambda x, y: u(x, y, t)
        v_t = lambda x, y: v(x, y, t)
        image = MyImage.image_from_matrix(video[x_min:x_max+1, y_min:y_max+1, t], root + "motion_stop" + "_" + str(t) + ".bmp")
        image.shift_image(u_t, v_t, root + "motion_stop" + "_" + str(t) + ".bmp")


def wave(args):
    if args["w"] is None:
        return

    video = args["video"]
    frequency, pixel_size, x_min, y_min, x_max, y_max = args["w"]
    if x_max <= 0:
        x_max += video.width-1
    if y_max <= 0:
        y_max += video.height-1
    
    quality_level = args["q"]
    data = args["data"]

    wave_data = data[int(x_min): int(x_max)+1, int(y_min)+1: int(y_max)+1, :]

    args["wave_data"] = wave_data
    args["frequency"] = frequency
    args["pixel_size"] = pixel_size


DESCRIPTION = "Given a sequence of 8 images, this program applies a combined method of the algorithms presented in the Timoner-Freeman and Horn-Schnuck papers in order to " + \
              "compute the optical flow at each pixel and from frame to frame. Then, these 7 results are fit into a sine wave in order to get the amplitude and phase of the motion at each pixel.\n" + \
              "For computing the optical flow, the first step consists in computing the gradients with respect to the x-direction, y-direction and time with the method described in the Timoner-Freeman paper. " + \
              "The conventions are that the x-axis is horizontal and the y-axis is vertical to the image, where the point (x, y) = (0, 0) represents the most upper-left pixel. " + \
              "The second step consists in dividing the image into squares of size win_min x win_min. For each of them, we consider the gradients of all the pixels in a square of a bigger size " + \
              "determined by the parameter win_max. These gradients are used to determine the displacement of this square from one frame to another. " + \
              "For these calculations, this algorithm uses a 2x2 matrix whose eigenvalues are good for corner dectection. Corners, defined as regions of good contrast in both " + \
              "x-direction and y-direction, contain reliable answers. Therefore, this algorithm uses a parameter called quality-level (a float from 0 to 1) " + \
              "which ignores all the squares whose ratio between the biggest two eigenvalues and its own is smaller than or equal to this parameter. " + \
              "After determining which squares are good, the centers of the good squares are initialized with the computed displacements and then used to linearly interpolate to get the displacements " + \
              "at each pixel of the image from one frame to the next.\n" +\
              "The third and final step is to use the displacements computed in the second step as the initial values for the Horn-Schunk algorithm. This algorithm depends on two variables. The first one is " +\
              "the smoothness parameter. The higher this value is, the smoother but less accurate the answers would be. However, if the smoothness parameter is too low, then the algorithm is very noise sensitive." + \
              "The second parameter is the number of iterations. Due to the second step, less iterations are needed since the Timoner-Freeman algorithm is very accurate in regions with good gradients. " + \
              "In addition, this program modifies the original Horn-Schunk algorithm by using the eigenvalues previously computed in the second step as weights. The higher the weight, the more similar the final value is " + \
              "to the displacements computed in the second step.\n" +\
              "The default behavior of this program consists in using the parameters win_min=5, win_max=25, quality_level=0.07, max_iterations=100, smoothness=100 to compute the optical flow and then " + \
              "use this information to calculate the amplitudes and phases in the x and y direction after fitting them to a sinousoidal movement. This programs generates four heat maps plots containing this information as well as " + \
              "saving them as png files, whose names are #amplitudes_x.png, #amplitudes_y.png, #phases_x.png, #phases_y.png, where # is the positional argument root." 


HELP = {"image": "8 file paths of the images to be analysed. After sorting them by their name, the i-th image will be consider as the i-th frame in the motion analysis",
        "root": "the root used to store all the files to be saved",
        "-c": "specify a region of interest, where x_min and y_min represent the top-left corner of the ROI and x_max and y_max " + \
              "represent the bottom right corner of the ROI. For x_max and y_max, non-positive integers are allowed and, in that case, the value " + \
              "resulting from substracting from the width-1 and the height-1 of the image will be used instead",
        "-x": "save the csv files #phases_x.csv, #phases_y.csv, #amplitudes_x.csv, #amplitudes_y.csv, where # is the positional argument root, containing the " + \
              "phases in the x-direction, y-direction and amplitudes in the x-direction and y-direction of every pixel in " + \
              "the given ROI after fitting the motion to a sine wave. The cell at column x and row y represents the pixel" + \
              "(x, y) relative to the ROI",
        "-s": "disable option of saving all plots as png's",
        "-p": "disable option of showing the plots at the end of the execution of the program",
        "-t": "only display amplitudes and phases in x and y of pixels with amplitudes in x and y greater than thr_x and thr_y",
        "-q": "quality-level parameter of this algorithm as described above",
        "-w": "outputs wave speed and decay constant given a ROI. The parameter frequency is in HZ and the size of a pixel in nanometers. x_min, y_min, x_max, y_max work " + \
              "analogously as for flag -c, but they are relative to the ROI specified with -c if used. The decay constant is in 1/nm and the wave speed in m/s.",
        "-a": "plots the amplitudes as a vector field instead as heat maps. By using this flag, the phases and amplitudes heat maps will neither be shown nor saved. " + \
              "The parameter k is an integer that indicates the distance between each arrow of the vector field. The parameter scale is a positive float which determines by how much "
              "the lenghts of the arrows will be multiplied. For instance, if scale is 10, the arrows in the vector field will be 10 bigger than their original value. " + \
              "This plot just appears at the end of the execution of this program; no file is saved.",
        "--a8": "save 7 png files in which the i-th file (i from 0 to 6) consists of a vector field representing the optical flow from frame i to frame i+1. "+\
                "The i-th file is named #displacements_vector_field_i.png, where # is the positional argument root. The parameters k and scale work as described in flag -a",
        "--hornShunck": "run Horn and Shunck algorithm if specified",
        "--fullImage": "run algoritm in full image. If not set, downsample image",
        "--includeDisplacements": "outputs frame to frame displacements",
        "--includeErrors": "outputs sine fitting error at each pixel",
        "--motionmag": "store 8 images resulting from the motion magnification of the original 8 frames. The input factor determines " + \
                       "the factor by which the displacements will be multiplied. The name of the file of the i-th frame (with i between 0 and 7) will be #motion_mag_factorF_i.bmp, where # is the " + \
                       "positional argument root and F is the float of the input factor.",
        "--motionstop": "shift all the images by negating the displacements computed. If the algorithm is accurate, " + \
                        "all the images should be similar to the first frame. The name of the file of the i-th frame (with i between 0 and 7) " +\
                        "will be #motion_stop_i.bmp"
        }

if __name__ == "__main__":
    flags = [{"image" : dict(action="store", nargs=8, help=HELP["image"], metavar="frame")},
             {"root": dict(action="store", help=HELP["root"])},
             {"-c": dict(action="store", help=HELP["-c"], nargs=4, metavar=("x_min", "y_min", "x_max", "y_max"), type=int, default=(0,0,0,0))},
             {"-q": dict(action="store", help=HELP["-q"], metavar="quality-level", type=float, default=0.07)},
             {"-t": dict(action="store", help=HELP["-t"], nargs=2, metavar=("thr_x", "thr_y"), default=(0.0, 0.0), type=float)},
             {"-a": dict(action="store", help=HELP["-a"], nargs=2, metavar=("k", "scale"), type=float)},
             {"--a8": dict(action="store", help=HELP["--a8"], nargs=2, metavar=("k", "scale"), type=float)},
             {"-s": dict(action="store_false", help=HELP["-s"])},
             {"-p": dict(action="store_false", help=HELP["-p"])},
             {"-x": dict(action="store_true", help=HELP["-x"])},
             {"-w": dict(action="store", help=HELP["-w"], nargs=6, metavar=("frequency", "pixel_dimensions", "x_min", "y_min", "x_max", "y_max"), type=float)},
             {"--hornShunck": dict(action="store_true", help=HELP["--hornShunck"])},
             {"--fullImage": dict(action="store_true", help=HELP["--fullImage"])},
             {"--includeDisplacements": dict(action="store_true", help=HELP["--includeDisplacements"])},
             {"--includeErrors": dict(action="store_true", help=HELP["--includeErrors"])},
             {"--motionmag": dict(action="store", help=HELP["--motionmag"], type=float, metavar="factor")},
             {"--motionstop": dict(action="store_true", help=HELP["--motionstop"])}]

    functions = [generate_data, save_data, motion_mag, motion_stop, wave, plot_data]

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    for flag in flags:
        key = list(flag.keys())[0]
        parser.add_argument(key, **flag[key])

    args = vars(parser.parse_args(sys.argv[1:]))

    for f in functions:
        f(args)




