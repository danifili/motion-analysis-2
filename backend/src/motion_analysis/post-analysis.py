import argparse
import numpy as np
import sys
import matplotlib.pyplot as plt
from src.motion_analysis.Plot import Plot
from src.motion_analysis.MyImage import MyImage


def generate_data(args):
    amplitudes_x = np.genfromtxt(args['amplitudes_x'], delimiter=',')
    amplitudes_y = np.genfromtxt(args['amplitudes_y'], delimiter=',')
    phases_x = np.genfromtxt(args['phases_x'], delimiter=',')
    phases_y = np.genfromtxt(args['phases_y'], delimiter=',')
    data = np.zeros(amplitudes_x.shape + (4,))
    for i, matrix in enumerate([amplitudes_x, amplitudes_y, phases_x, phases_y]):
        data[:, :, i] = matrix
    return data


def apply_t(data, args):
    thr_x, thr_y = args['t']
    mask_x = data[:, :, 0] < thr_x
    mask_y = data[:, :, 1] < thr_y
    data[:, :, 0][mask_x] = 'nan'
    data[:, :, 2][mask_x] = 'nan'
    data[:, :, 1][mask_y] = 'nan'
    data[:, :, 3][mask_y] = 'nan'


def show_plot(data, args, image):
    min_corner = args['c'][:2]
    max_corner = args['c'][2:4]
    max_corner = data.shape[:2][0] - max_corner[0], data.shape[:2][1] - max_corner[1]
    x_min, y_min = min_corner
    x_max, y_max = max_corner
    matrix = np.zeros((data.shape[:2]) + (1,))
    matrix[:, :, 0] = image[0:data.shape[:2][0], 0:data.shape[:2][1]]
    max_corner = x_max-1, y_max-1

    amplitudes_x_fig = plt.figure('Amplitudes X')
    plt.title('Amplitudes X')
    Plot.scalar_heat_map(matrix, min_corner, max_corner, 0, data[x_min:x_max,y_min:y_max,0], alpha=0.3)

    amplitudes_y_fig = plt.figure('Amplitudes Y')
    plt.title('Amplitudes Y')
    Plot.scalar_heat_map(matrix, min_corner, max_corner, 0, data[x_min:x_max,y_min:y_max,1], alpha=0.3)

    phases_x_fig = plt.figure('Phases X')
    plt.title('Phases X')
    Plot.phase_heat_map(matrix, min_corner, max_corner, 0, data[x_min:x_max,y_min:y_max,2], alpha=0.3)

    phases_y_fig = plt.figure('Phases Y')
    plt.title('Phases Y')
    Plot.phase_heat_map(matrix, min_corner, max_corner, 0, data[x_min:x_max,y_min:y_max,3], alpha=0.3)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = [{'image': dict(action="store",)},
             {"amplitudes_x": dict(action="store",)},
             {"amplitudes_y": dict(action="store",)},
             {"phases_x": dict(action="store",)},
             {"phases_y": dict(action="store",)},
             {"-s": dict(action="store_false",)},
             {"-p": dict(action="store_false",)},
             {"-t": dict(action="store",nargs=2, metavar=("thr_x", "thr_y"), default=(0.0, 0.0),type=float)},
             {"-c": dict(action="store", nargs=4, metavar=("x_min", "y_min", "x_max", "y_max"), type=int, default=(0, 0, 0, 0))}]

    modifiers_callbacks = [("t", apply_t)]
    action_callbacks = [("t", show_plot)]

    for flag in flags:
        key = list(flag.keys())[0]
        parser.add_argument(key, **flag[key])

    args = vars(parser.parse_args(sys.argv[1:]))

    data = generate_data(args)
    print (args['image'])
    image = MyImage(args['image'])

    for flag, callback in modifiers_callbacks:
        if args[flag]:
            callback(data, args)

    for flag, callback in action_callbacks:
        if args[flag]:
            callback(data, args, image)


