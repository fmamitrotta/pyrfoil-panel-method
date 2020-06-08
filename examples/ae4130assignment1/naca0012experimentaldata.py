"""Experimental data of NACA 0012

This module allows the user to retrieve the experimental data of NACA 0012. Such experimental data consists of both
force and pressure data. Everything is taken from https://turbmodels.larc.nasa.gov/naca0012_val.html.

This script requires that `numpy` be installed within the Python environment you are running this script in.

This file contains the following functions:

    * abbott_force_data - returns data from Abbott & Von Doenhoff, "Theory of Wing Sections"
    * ladson_force_data - returns data from Ladson, NASA TM 4074, 1988
    * ladson_pressure_data - returns data from Ladson, Hill, & Johnson, NASA TM 100526, 1987
    * gregory_pressure_data - returns data from Gregory & O'Reilly, NASA R&M 3726, Jan 1970
"""

import numpy as np


def abbott_force_data():
    # Digitized data from Abbott & Von Doenhoff,
    # "Theory of Wing Sections"
    # NACA 0012 wing section, Re=6 million
    # (digitizing only approximate, due to poor quality
    # of plot in the book)
    # variables = "alpha, deg", "cl"
    return np.array([[-17.2794, -1.25323],
                     [-16.2296, -1.34704],
                     [-15.8616, -1.54416],
                     [-15.1713, -1.51805],
                     [-14.3133, -1.44038],
                     [-13.2811, -1.3712],
                     [-12.2535, -1.25912],
                     [-11.2222, -1.18135],
                     [-10.1947, -1.06927],
                     [-8.14138, -0.827958],
                     [-6.25579, -0.638207],
                     [-5.22822, -0.526128],
                     [-4.19972, -0.422627],
                     [-1.96944, -0.215533],
                     [0., 0.],
                     [0.940006, 0.120611],
                     [1.96944, 0.215533],
                     [2.99515, 0.34477],
                     [3.85131, 0.439599],
                     [4.87888, 0.551678],
                     [5.90831, 0.6466],
                     [7.96346, 0.870758],
                     [10.1891, 1.12074],
                     [11.0471, 1.19842],
                     [13.1088, 1.36252],
                     [16.3759, 1.59591],
                     [16.5678, 1.42443],
                     [17.2971, 1.09024]])


def ladson_force_data():
    # Data from Ladson, NASA TM 4074, 1988
    # Re=6 million, transition fixed with different  no. of grit
    # M=0.15
    # variables = "alpha, deg", "cl", "cd"
    return {
        '80 grit': np.array([[-4.04, -.4417, .00871],
                             [-2.14, -.2385, .00800],
                             [-.05, -.0126, .00809],
                             [2.05, .2125, .00816],
                             [4.04, .4316, .00823],
                             [6.09, .6546, .00885],
                             [8.30, .8873, .01050],
                             [10.12, 1.0707, .01201],
                             [11.13, 1.1685, .01239],
                             [12.12, 1.2605, .01332],
                             [13.08, 1.3455, .01503],
                             [14.22, 1.4365, .01625],
                             [15.26, 1.5129, .01900],
                             [16.30, 1.5739, .02218],
                             [17.13, 1.6116, .02560],
                             [18.02, .9967, .18785],
                             [19.08, 1.1358, .27292]]),
        '120 grit': np.array([[-4.01, -.4466, .00843],
                              [-2.12, -.2425, .00789],
                              [-.01, -.0120, .00811],
                              [.01, -.0122, .00804],
                              [2.15, .2236, .00823],
                              [4.11, .4397, .00879],
                              [6.01, .6487, .00842],
                              [8.08, .8701, .00995],
                              [10.10, 1.0775, .01175],
                              [11.23, 1.1849, .01248],
                              [12.13, 1.2720, .01282],
                              [13.26, 1.3699, .01408],
                              [14.30, 1.4571, .01628],
                              [15.27, 1.5280, .01790],
                              [16.16, 1.5838, .02093],
                              [17.24, 1.6347, .02519],
                              [18.18, 1.1886, .25194],
                              [19.25, 1.1888, .28015]]),
        '180 grit': np.array([[-3.99, -.4363, .00871],
                              [-1.98, -.2213, .00792],
                              [-.03, -.0115, .00803],
                              [.04, -.0013, .00811],
                              [2.00, .2213, .00814],
                              [4.06, .4365, .00814],
                              [6.09, .6558, .00851],
                              [8.09, .8689, .00985],
                              [10.18, 1.0809, .01165],
                              [11.13, 1.1731, .01247],
                              [12.10, 1.2644, .01299],
                              [13.31, 1.3676, .01408],
                              [14.08, 1.4316, .01533],
                              [15.24, 1.5169, .01870],
                              [16.33, 1.5855, .02186],
                              [17.13, 1.6219, .02513],
                              [18.21, 1.0104, .25899],
                              [19.27, 1.0664, .43446]])
    }


def ladson_pressure_data():
    # Data from Ladson, Hill, & Johnson, NASA TM 100526, 1987
    # Re=3 million
    # M=0.30
    # Note: may not be sufficiently 2D, because aspect ratio of model only 1.333
    # variables = "x/c", "cp"
    return {
        'alpha=0': np.array([[.9483, .0657],
                             [.9000, -.0013],
                             [.8503, -.0383],
                             [.7998, -.0679],
                             [.7497, -.1128],
                             [.7003, -.1324],
                             [.6502, -.1682],
                             [.5997, -.1978],
                             [.5506, -.2222],
                             [.5000, -.2520],
                             [.4503, -.2730],
                             [.4000, -.3098],
                             [.3507, -.3320],
                             [.3002, -.3603],
                             [.2501, -.3904],
                             [.2004, -.4151],
                             [.1504, -.4400],
                             [.1000, -.4375],
                             [.0755, -.4112],
                             [.0510, -.3830],
                             [.0251, -.2626],
                             [.0122, -.0181],
                             [.0000, 1.0029],
                             [.0135, -.0047],
                             [.0271, -.2347],
                             [.0515, -.3617],
                             [.0763, -.4026],
                             [.1012, -.4263],
                             [.1503, -.4106],
                             [.1994, -.4022],
                             [.2501, -.3760],
                             [.2999, -.3516],
                             [.3499, -.3201],
                             [.3994, -.3007],
                             [.4496, -.2668],
                             [.4997, -.2449],
                             [.5492, -.2077],
                             [.5994, -.1885],
                             [.6495, -.1547],
                             [.6996, -.1279],
                             [.7489, -.1093],
                             [.8003, -.0659],
                             [.8500, -.0347],
                             [.8993, .0158],
                             [.9489, .0597]]),
        'alpha=10': np.array([[.9483, .0691],
                              [.9000, .0335],
                              [.8503, .0402],
                              [.7998, .0416],
                              [.7497, .0303],
                              [.7003, .0389],
                              [.6502, .0389],
                              [.5997, .0413],
                              [.5506, .0397],
                              [.5000, .0469],
                              [.4503, .0627],
                              [.4000, .0751],
                              [.3507, .0989],
                              [.3002, .1311],
                              [.2501, .1686],
                              [.2004, .2294],
                              [.1504, .3276],
                              [.1000, .4774],
                              [.0755, .5914],
                              [.0510, .7384],
                              [.0251, .9543],
                              [.0122, .9863],
                              [.0000, -2.3938],
                              [.0135, -3.8803],
                              [.0271, -2.9668],
                              [.0515, -2.3162],
                              [.0763, -1.9884],
                              [.1012, -1.7245],
                              [.1503, -1.3903],
                              [.1994, -1.1922],
                              [.2501, -1.0278],
                              [.2999, -.8996],
                              [.3499, -.7865],
                              [.3994, -.6937],
                              [.4496, -.6044],
                              [.4997, -.5364],
                              [.5492, -.4546],
                              [.5994, -.3846],
                              [.6495, -.3210],
                              [.6996, -.2626],
                              [.7489, -.2153],
                              [.8003, -.1460],
                              [.8500, -.0852],
                              [.8993, -.0178],
                              [.9489, .0475]]),
        'alpha=15': np.array([[.9483, -.0641],
                              [.9000, -0.573],
                              [.8503, -.0203],
                              [.7998, .0098],
                              [.7497, .0216],
                              [.7003, .0398],
                              [.6502, .0539],
                              [.5997, .0813],
                              [.5506, .0936],
                              [.5000, .1240],
                              [.4503, .1591],
                              [.4000, .1797],
                              [.3507, .2172],
                              [.3002, .2785],
                              [.2501, .3383],
                              [.2004, .4161],
                              [.1504, .5226],
                              [.1000, .6968],
                              [.0755, .7988],
                              [.0510, .9190],
                              [.0251, .9931],
                              [.0122, .7638],
                              [.0000, -5.1051],
                              [.0135, -5.9032],
                              [.0271, -4.0778],
                              [.0515, -2.9754],
                              [.0763, -2.4752],
                              [.1012, -2.1222],
                              [.1503, -1.6276],
                              [.1994, -1.3478],
                              [.2501, -1.1165],
                              [.2999, -.9549],
                              [.3499, -.8456],
                              [.3994, -.7484],
                              [.4496, -.6671],
                              [.4997, -.5969],
                              [.5492, -.5463],
                              [.5994, -.4974],
                              [.6495, -.4553],
                              [.6996, -.4128],
                              [.7489, -.3733],
                              [.8003, -.3226],
                              [.8500, -.2794],
                              [.8993, .3055],
                              [.9489, -.1914]]),
        'alpha=0, fixed transition': np.array([[.9463, .1336],
                                               [.9000, .0631],
                                               [.8503, .232],
                                               [.7998, -.0126],
                                               [.7497, -.0388],
                                               [.7003, -.0762],
                                               [.6502, -.1036],
                                               [.5997, -.1265],
                                               [.5506, -.1573],
                                               [.5000, -.1851],
                                               [.4503, -.2094],
                                               [.4000, -.2355],
                                               [.3507, -.2692],
                                               [.3002, -.2949],
                                               [.2501, -.3258],
                                               [.2004, -.3573],
                                               [.1504, -.3621],
                                               [.1000, -.3604],
                                               [.0755, -.3427],
                                               [.0510, -.2886],
                                               [.0251, -.1440],
                                               [.0122, .0721],
                                               [.0000, 1.0561],
                                               [.0135, .0468],
                                               [.0271, -.1666],
                                               [.0515, -.2812],
                                               [.0763, -.3477],
                                               [.1012, -.3599],
                                               [.1503, -.3438],
                                               [.1994, -.3528],
                                               [.2501, -.3210],
                                               [.2999, -.2962],
                                               [.3499, -.2677],
                                               [.3994, -.2345],
                                               [.4496, -.2104],
                                               [.4997, -.1824],
                                               [.5492, -.1493],
                                               [.5994, -.1239],
                                               [.6495, -.0952],
                                               [.6996, -.0755],
                                               [.7489, -.0409],
                                               [.8003, -.0084],
                                               [.8500, .0263],
                                               [.8993, .0740],
                                               [.9489, .1285]]),
        'alpha=10, fixed transition': np.array([[.9483, .1335],
                                                [.9000, .1026],
                                                [.8503, .1037],
                                                [.7998, .1057],
                                                [.7497, .1031],
                                                [.7003, .1056],
                                                [.6502, .0920],
                                                [.5997, .1147],
                                                [.5506, .1132],
                                                [.5000, .1261],
                                                [.4503, .1424],
                                                [.4000, .1597],
                                                [.3507, .1786],
                                                [.3002, .2145],
                                                [.2501, .2537],
                                                [.2004, .3154],
                                                [.1504, .4047],
                                                [.1000, .5591],
                                                [.0755, .6593],
                                                [.0510, .8124],
                                                [.0251, 1.0204],
                                                [.0122, 1.0505],
                                                [.0000, -2.3905],
                                                [.0135, -3.8320],
                                                [.0271, -2.8263],
                                                [.0515, -2.3446],
                                                [.0763, -1.9334],
                                                [.1012, -1.6899],
                                                [.1503, -1.3465],
                                                [.1994, -1.1389],
                                                [.2501, -.9749],
                                                [.2999, -.8345],
                                                [.3499, -.7279],
                                                [.3994, -.6274],
                                                [.4496, -.5438],
                                                [.4997, -.4624],
                                                [.5492, -.3950],
                                                [.5994, -.3236],
                                                [.6495, -.2723],
                                                [.6996, -.2003],
                                                [.7489, -.1446],
                                                [.8003, -.0766],
                                                [.8500, -.0176],
                                                [.8993, .0468],
                                                [.9489, .1108]]),
        'alpha=15, fixed transition': np.array([[.9483, .0769],
                                                [.9000, .0591],
                                                [.8503, .0904],
                                                [.7998, .1027],
                                                [.7497, .1199],
                                                [.7003, .1397],
                                                [.6502, .1612],
                                                [.5997, .1800],
                                                [.5506, .1913],
                                                [.5000, .2270],
                                                [.4503, .2372],
                                                [.4000, .2774],
                                                [.3507, .3155],
                                                [.3002, .3586],
                                                [.2501, .4368],
                                                [.2004, .5185],
                                                [.1504, .6214],
                                                [.1000, .7971],
                                                [.0755, .8800],
                                                [.0510, 1.0054],
                                                [.0251, 1.0296],
                                                [.0122, .7393],
                                                [0.000, -5.8153],
                                                [.0135, -6.4597],
                                                [.0271, -4.2582],
                                                [.0515, -3.2200],
                                                [.0763, -2.6220],
                                                [.1012, -2.2710],
                                                [.1503, -1.7264],
                                                [.1994, -1.4306],
                                                [.2501, -1.1843],
                                                [.2999, -.9787],
                                                [.3499, -.8619],
                                                [.3994, -.7367],
                                                [.4496, -.6374],
                                                [.4997, -.5194],
                                                [.5492, -.4446],
                                                [.5994, -.3623],
                                                [.6495, -.2907],
                                                [.6996, -.2363],
                                                [.7489, -.1862],
                                                [.8003, -.1372],
                                                [.8500, -.0831],
                                                [.8993, -.0601],
                                                [.9489, -.0016]])
    }


def gregory_pressure_data():
    # Data from Gregory & O'Reilly, NASA R&M 3726, Jan 1970
    # Re=2.88 million, free transition
    # Data was digitized from a photocopy - hence is only approximate
    # Data for upper airfoil surface only
    # variables = "x/c", "cp"
    return {
        'alpha=0': np.array([[0, 1],
                             [0.0023497, 0.847673],
                             [0.00496048, 0.456198],
                             [0.00526903, 0.173569],
                             [0.0142406, -0.044407],
                             [0.0209337, -0.175278],
                             [0.0473501, -0.372653],
                             [0.0779437, -0.396388],
                             [0.0976194, -0.41941],
                             [0.128166, -0.418874],
                             [0.150001, -0.411087],
                             [0.178387, -0.402938],
                             [0.289702, -0.36672],
                             [0.322431, -0.347115],
                             [0.387891, -0.307906],
                             [0.448983, -0.268412],
                             [0.514442, -0.229203],
                             [0.579902, -0.189994],
                             [0.638834, -0.159098],
                             [0.704317, -0.114629],
                             [0.767593, -0.065278],
                             [0.835236, -0.026211],
                             [0.896305, 0.03502],
                             [0.959533, 0.0978565],
                             [1.0009, 0.173854]]),
        'alpha=10': np.array([[0, -3.66423],
                              [0.00218341, -5.04375],
                              [0.00873362, -5.24068],
                              [0.0131004, -4.67125],
                              [0.0174672, -4.32079],
                              [0.0480349, -2.74347],
                              [0.0742358, -2.26115],
                              [0.0982533, -1.95405],
                              [0.124454, -1.7345],
                              [0.146288, -1.55884],
                              [0.176856, -1.36109],
                              [0.28821, -1.00829],
                              [0.320961, -0.941877],
                              [0.384279, -0.787206],
                              [0.447598, -0.654432],
                              [0.515284, -0.543461],
                              [0.576419, -0.432633],
                              [0.637555, -0.343703],
                              [0.700873, -0.254725],
                              [0.766376, -0.1657],
                              [0.831878, -0.098572],
                              [0.893013, -0.00964205],
                              [0.958515, 0.0793835],
                              [1, 0.124088]]),
        'alpha=15': np.array([[-7.59438e-05, -8.65066],
                              [0.0024302, -10.1789],
                              [0.00450442, -9.72033],
                              [0.00870506, -9.04329],
                              [0.0129722, -8.67192],
                              [0.0167741, -6.16084],
                              [0.0467387, -3.99796],
                              [0.0769928, -3.16694],
                              [0.0964534, -2.68574],
                              [0.146315, -2.05038],
                              [0.174528, -1.83081],
                              [0.287443, -1.23636],
                              [0.317853, -1.12586],
                              [0.380854, -0.9266],
                              [0.443854, -0.727343],
                              [0.509042, -0.593492],
                              [0.576404, -0.459546],
                              [0.635076, -0.347813],
                              [0.698095, -0.235891],
                              [0.761123, -0.167637],
                              [0.8285, -0.0991921],
                              [0.893707, -0.0526765],
                              [0.954576, -0.0500185],
                              [1.00022, -0.00435728]])
    }
