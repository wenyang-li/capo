import numpy as np
east_x = 0.0050833611111080756
south_x = -6.4968000628571421
east_y = -3.1208064055555544
south_y = 1.8538589797142855
antpos = {
    'nant': 128,
    #EastHex
        57: {'top_x':  86.352 + east_x,  'top_y':  287.221067 + east_y},
        58: {'top_x': 100.352 + east_x,  'top_y':  287.221067 + east_y},
        59: {'top_x': 114.352 + east_x,  'top_y':  287.221067 + east_y},
        60: {'top_x': 128.352 + east_x,  'top_y':  287.221067 + east_y},
        
        61: {'top_x':  79.352 + east_x,  'top_y':  275.0967113 + east_y},
        62: {'top_x':  93.352 + east_x,  'top_y':  275.0967113 + east_y},
        63: {'top_x': 107.352 + east_x,  'top_y':  275.0967113 + east_y},
        64: {'top_x': 121.352 + east_x,  'top_y':  275.0967113 + east_y},
        65: {'top_x': 135.352 + east_x,  'top_y':  275.0967113 + east_y},
        
        66: {'top_x':  72.352 + east_x,  'top_y':  262.9723557 + east_y},
        67: {'top_x':  86.352 + east_x,  'top_y':  262.9723557 + east_y},
        68: {'top_x': 100.352 + east_x,  'top_y':  262.9723557 + east_y},
        69: {'top_x': 114.352 + east_x,  'top_y':  262.9723557 + east_y},
        70: {'top_x': 128.352 + east_x,  'top_y':  262.9723557 + east_y},
        71: {'top_x': 142.352 + east_x,  'top_y':  262.9723557 + east_y},
        
        72: {'top_x':  65.352 + east_x,  'top_y': 250.848 + east_y},
        73: {'top_x':  79.352 + east_x,  'top_y': 250.848 + east_y},
        74: {'top_x':  93.352 + east_x,  'top_y': 250.848 + east_y},
        75: {'top_x': 121.352 + east_x,  'top_y': 250.848 + east_y},
        76: {'top_x': 135.352 + east_x,  'top_y': 250.848 + east_y},
        77: {'top_x': 149.352 + east_x,  'top_y': 250.848 + east_y},
        
        78: {'top_x':  72.352 + east_x,  'top_y': 238.7236444 + east_y},
        79: {'top_x':  86.352 + east_x,  'top_y': 238.7236444 + east_y},
        80: {'top_x': 100.352 + east_x,  'top_y': 238.7236444 + east_y},
        81: {'top_x': 114.352 + east_x,  'top_y': 238.7236444 + east_y},
        82: {'top_x': 128.352 + east_x,  'top_y': 238.7236444 + east_y},
        83: {'top_x': 142.352 + east_x,  'top_y': 238.7236444 + east_y},
        
        84: {'top_x':  79.352 + east_x,  'top_y': 226.5992887 + east_y},
        85: {'top_x':  93.352 + east_x,  'top_y': 226.5992887 + east_y},
        86: {'top_x': 107.352 + east_x,  'top_y': 226.5992887 + east_y},
        87: {'top_x': 121.352 + east_x,  'top_y': 226.5992887 + east_y},
        88: {'top_x': 135.352 + east_x,  'top_y': 226.5992887 + east_y},
        
        89: {'top_x':  86.352 + east_x,  'top_y': 214.474933 + east_y},
        90: {'top_x': 100.352 + east_x,  'top_y': 214.474933 + east_y},
        91: {'top_x': 114.352 + east_x,  'top_y': 214.474933 + east_y},
        92: {'top_x': 128.352 + east_x,  'top_y': 214.474933 + east_y},
        
        #south hex
        
#        92: {'top_x':   -3.49,  'top_y': 156.590067,  'cable':    90},
        93: {'top_x':   10.51 + south_x,  'top_y': 156.590067 + south_y},
        94: {'top_x':   24.51 + south_x,  'top_y': 156.590067 + south_y},
        95: {'top_x':   38.51 + south_x,  'top_y': 156.590067 + south_y},
        
        96: {'top_x':  -10.49 + south_x,  'top_y': 144.4657113 + south_y},
        97: {'top_x':    3.51 + south_x,  'top_y': 144.4657113 + south_y},
        98: {'top_x':   17.51 + south_x,  'top_y': 144.4657113 + south_y},
        99: {'top_x':   31.51 + south_x,  'top_y': 144.4657113 + south_y},
        100: {'top_x':   45.51 + south_x,  'top_y': 144.4657113 + south_y},
        
        101: {'top_x':  -17.49 + south_x,  'top_y': 132.3413557 + south_y},
        102: {'top_x':   -3.49 + south_x,  'top_y': 132.3413557 + south_y},
        103: {'top_x':   10.51 + south_x,  'top_y': 132.3413557 + south_y},
        104: {'top_x':   24.51 + south_x,  'top_y': 132.3413557 + south_y},
        105: {'top_x':   38.51 + south_x,  'top_y': 132.3413557 + south_y},
        106: {'top_x':   52.51 + south_x,  'top_y': 132.3413557 + south_y},
        
        107: {'top_x':  -24.49 + south_x,  'top_y': 120.217 + south_y},
        108: {'top_x':  -10.49 + south_x,  'top_y': 120.217 + south_y},
        109: {'top_x':    3.51 + south_x,  'top_y': 120.217 + south_y},
        110: {'top_x':   31.51 + south_x,  'top_y': 120.217 + south_y},
        111: {'top_x':   45.51 + south_x,  'top_y': 120.217 + south_y},
        112: {'top_x':   59.51 + south_x,  'top_y': 120.217 + south_y},
        
        113: {'top_x':  -17.49 + south_x,  'top_y': 108.0926444 + south_y},
        114: {'top_x':   -3.49 + south_x,  'top_y': 108.0926444 + south_y},
        115: {'top_x':   10.51 + south_x,  'top_y': 108.0926444 + south_y},
        116: {'top_x':   24.51 + south_x,  'top_y': 108.0926444 + south_y},
        117: {'top_x':   38.51 + south_x,  'top_y': 108.0926444 + south_y},
        118: {'top_x':   52.51 + south_x,  'top_y': 108.0926444 + south_y},
        
        119: {'top_x':  -10.49 + south_x,  'top_y': 95.96828869 + south_y},
        120: {'top_x':    3.51 + south_x,  'top_y': 95.96828869 + south_y},
        121: {'top_x':   17.51 + south_x,  'top_y': 95.96828869 + south_y},
        122: {'top_x':   31.51 + south_x,  'top_y': 95.96828869 + south_y},
        123: {'top_x':   45.51 + south_x,  'top_y': 95.96828869 + south_y},
        
        124: {'top_x':   -3.49 + south_x,  'top_y': 83.84393304 + south_y},
        125: {'top_x':   10.51 + south_x,  'top_y': 83.84393304 + south_y},
        126: {'top_x':   24.51 + south_x,  'top_y': 83.84393304 + south_y},
        127: {'top_x':   38.51 + south_x,  'top_y': 83.84393304 + south_y},
}

realpos={
 0: {'top_x': -149.785, 'top_y': 265.814, 'top_z': 377.01099},
 1: {'top_x': -95.364998, 'top_y': 270.18799, 'top_z': 377.285},
 2: {'top_x': -88.512001, 'top_y': 266.021, 'top_z': 377.30899},
 3: {'top_x': -78.697998, 'top_y': 258.431, 'top_z': 377.29099},
 4: {'top_x': -86.940002, 'top_y': 248.452, 'top_z': 377.33899},
 5: {'top_x': -102.325, 'top_y': 244.088, 'top_z': 377.258},
 6: {'top_x': -93.346001, 'top_y': 242.27901, 'top_z': 377.26999},
 7: {'top_x': -148.33, 'top_y': 220.271, 'top_z': 377.21399},
 8: {'top_x': -132.939, 'top_y': 282.61099, 'top_z': 377.00201},
 9: {'top_x': -139.867, 'top_y': 275.49799, 'top_z': 377.01801},
 10: {'top_x': -123.667, 'top_y': 284.854, 'top_z': 377.04999},
 11: {'top_x': -92.889, 'top_y': 281.358, 'top_z': 377.254},
 12: {'top_x': -77.834999, 'top_y': 282.83899, 'top_z': 377.25299},
 13: {'top_x': -83.475998, 'top_y': 277.56299, 'top_z': 377.29999},
 14: {'top_x': -76.987999, 'top_y': 272.63699, 'top_z': 377.298},
 15: {'top_x': -176.006, 'top_y': 289.90302, 'top_z': 376.62299},
 16: {'top_x': -23.620001, 'top_y': 201.54601, 'top_z': 376.80499},
 17: {'top_x': -23.802, 'top_y': 222.353, 'top_z': 376.87201},
 18: {'top_x': -40.488998, 'top_y': 235.17, 'top_z': 376.961},
 19: {'top_x': -54.527, 'top_y': 237.51801, 'top_z': 377.12299},
 20: {'top_x': -64.850998, 'top_y': 250.813, 'top_z': 377.23901},
 21: {'top_x': -70.074997, 'top_y': 259.40201, 'top_z': 377.26599},
 22: {'top_x': -60.056999, 'top_y': 263.67599, 'top_z': 377.17999},
 23: {'top_x': -50.675999, 'top_y': 264.54599, 'top_z': 377.10001},
 24: {'top_x': -67.989998, 'top_y': 271.103, 'top_z': 377.23599},
 25: {'top_x': -60.057999, 'top_y': 271.905, 'top_z': 377.20401},
 26: {'top_x': -70.885002, 'top_y': 278.81799, 'top_z': 377.25201},
 27: {'top_x': -45.540001, 'top_y': 273.245, 'top_z': 377.06601},
 28: {'top_x': -53.155998, 'top_y': 276.26801, 'top_z': 377.09},
 29: {'top_x': -22.659, 'top_y': 263.694, 'top_z': 376.948},
 30: {'top_x': -55.097, 'top_y': 284.285, 'top_z': 377.181},
 31: {'top_x': -14.9, 'top_y': 273.87, 'top_z': 376.91},
 32: {'top_x': -105.281, 'top_y': 217.17999, 'top_z': 377.255},
 33: {'top_x': -98.529999, 'top_y': 230.162, 'top_z': 377.26001},
 34: {'top_x': -81.845001, 'top_y': 229.33701, 'top_z': 377.26999},
 35: {'top_x': -79.541, 'top_y': 238.12801, 'top_z': 377.26999},
 36: {'top_x': -75.120003, 'top_y': 247.002, 'top_z': 377.255},
 37: {'top_x': -71.016998, 'top_y': 235.929, 'top_z': 377.23901},
 38: {'top_x': -62.567001, 'top_y': 228.987, 'top_z': 377.20999},
 39: {'top_x': -50.691002, 'top_y': 221.88499, 'top_z': 377.078},
 40: {'top_x': -160.465, 'top_y': 573.20697, 'top_z': 374.77301},
 41: {'top_x': -128.04601, 'top_y': 350.49399, 'top_z': 376.57001},
 42: {'top_x': -69.771004, 'top_y': 294.03101, 'top_z': 377.19199},
 43: {'top_x': -78.731003, 'top_y': 297.259, 'top_z': 377.19601},
 44: {'top_x': -103.941, 'top_y': 300.914, 'top_z': 377.04999},
 45: {'top_x': -100.172, 'top_y': 288.896, 'top_z': 377.15701},
 46: {'top_x': -395.452, 'top_y': 371.15201, 'top_z': 374.51501},
 47: {'top_x': -263.547, 'top_y': 389.27399, 'top_z': 375.32199},
 48: {'top_x': 326.14099, 'top_y': 203.608, 'top_z': 373.94501},
 49: {'top_x': 173.425, 'top_y': 193.595, 'top_z': 375.03201},
 50: {'top_x': -2.5150001, 'top_y': 292.841, 'top_z': 376.892},
 51: {'top_x': -46.554001, 'top_y': 287.46899, 'top_z': 377.11499},
 52: {'top_x': -56.140999, 'top_y': 295.922, 'top_z': 377.15601},
 53: {'top_x': -33.705002, 'top_y': 351.048, 'top_z': 376.97601},
 54: {'top_x': 50.126999, 'top_y': 536.159, 'top_z': 376.181},
 55: {'top_x': 84.879997, 'top_y': 519.14099, 'top_z': 376.314},
 56: {'top_x': 143.44099, 'top_y': 27.773001, 'top_z': 374.78},
 57: {'top_x': 86.373001, 'top_y': 284.121, 'top_z': 376.25},
 58: {'top_x': 100.347, 'top_y': 284.09698, 'top_z': 376.17999},
 59: {'top_x': 114.346, 'top_y': 284.09201, 'top_z': 375.91},
 60: {'top_x': 128.34, 'top_y': 284.104, 'top_z': 375.76999},
 61: {'top_x': 79.369003, 'top_y': 271.98599, 'top_z': 376.23999},
 62: {'top_x': 93.338997, 'top_y': 271.97601, 'top_z': 376.07001},
 63: {'top_x': 107.365, 'top_y': 271.957, 'top_z': 375.95001},
 64: {'top_x': 121.355, 'top_y': 271.991, 'top_z': 375.82999},
 65: {'top_x': 135.367, 'top_y': 271.96201, 'top_z': 375.64999},
 66: {'top_x': 72.379997, 'top_y': 259.84799, 'top_z': 376.23999},
 67: {'top_x': 86.348999, 'top_y': 259.84399, 'top_z': 376.10001},
 68: {'top_x': 100.373, 'top_y': 259.84698, 'top_z': 375.95999},
 69: {'top_x': 114.357, 'top_y': 259.83899, 'top_z': 375.82001},
 70: {'top_x': 128.34599, 'top_y': 259.84799, 'top_z': 375.70001},
 71: {'top_x': 142.36, 'top_y': 259.849, 'top_z': 375.5},
 72: {'top_x': 65.362, 'top_y': 247.72301, 'top_z': 376.23001},
 73: {'top_x': 79.365997, 'top_y': 247.742, 'top_z': 376.07999},
 74: {'top_x': 93.355003, 'top_y': 247.743, 'top_z': 375.94},
 75: {'top_x': 121.369, 'top_y': 247.744, 'top_z': 375.72},
 76: {'top_x': 135.353, 'top_y': 247.72099, 'top_z': 375.59},
 77: {'top_x': 149.358, 'top_y': 247.731, 'top_z': 375.41},
 78: {'top_x': 72.362999, 'top_y': 235.61, 'top_z': 376.09},
 79: {'top_x': 86.344002, 'top_y': 235.606, 'top_z': 375.97},
 80: {'top_x': 100.367, 'top_y': 235.616, 'top_z': 375.85001},
 81: {'top_x': 114.335, 'top_y': 235.59801, 'top_z': 375.73999},
 82: {'top_x': 128.349, 'top_y': 235.61, 'top_z': 375.56},
 83: {'top_x': 142.34, 'top_y': 235.586, 'top_z': 375.44},
 84: {'top_x': 79.374001, 'top_y': 223.474, 'top_z': 375.98999},
 85: {'top_x': 93.343002, 'top_y': 223.48801, 'top_z': 375.87},
 86: {'top_x': 107.357, 'top_y': 223.46201, 'top_z': 375.76999},
 87: {'top_x': 121.372, 'top_y': 223.48801, 'top_z': 375.57001},
 88: {'top_x': 135.35201, 'top_y': 223.48599, 'top_z': 375.48001},
 89: {'top_x': 86.362, 'top_y': 211.341, 'top_z': 375.88},
 90: {'top_x': 100.364, 'top_y': 211.338, 'top_z': 375.75},
 91: {'top_x': 114.359, 'top_y': 211.353, 'top_z': 375.60001},
 92: {'top_x': 128.345, 'top_y': 211.358, 'top_z': 375.48001},
 93: {'top_x': 4.006, 'top_y': 158.45599, 'top_z': 376.41},
 94: {'top_x': 18.021999, 'top_y': 158.439, 'top_z': 376.26001},
 95: {'top_x': 32.013, 'top_y': 158.442, 'top_z': 376.14001},
 96: {'top_x': -17.013, 'top_y': 146.30901, 'top_z': 376.59},
 97: {'top_x': -2.9849999, 'top_y': 146.31599, 'top_z': 376.42001},
 98: {'top_x': 11.022, 'top_y': 146.311, 'top_z': 376.31},
 99: {'top_x': 25.002001, 'top_y': 146.32001, 'top_z': 376.14001},
 100: {'top_x': 39.014, 'top_y': 146.304, 'top_z': 376.04001},
 101: {'top_x': -23.966, 'top_y': 134.201, 'top_z': 376.57001},
 102: {'top_x': -9.9829998, 'top_y': 134.185, 'top_z': 376.42001},
 103: {'top_x': 4.0050001, 'top_y': 134.20599, 'top_z': 376.29999},
 104: {'top_x': 17.99, 'top_y': 134.18201, 'top_z': 376.14999},
 105: {'top_x': 32.008999, 'top_y': 134.18401, 'top_z': 376.04999},
 106: {'top_x': 46.015999, 'top_y': 134.22099, 'top_z': 375.92001},
 107: {'top_x': -30.990999, 'top_y': 122.075, 'top_z': 376.60001},
 108: {'top_x': -16.989, 'top_y': 122.07, 'top_z': 376.45999},
 109: {'top_x': -2.9719999, 'top_y': 122.058, 'top_z': 376.32001},
 110: {'top_x': 25.014999, 'top_y': 122.055, 'top_z': 376.09},
 111: {'top_x': 39.014, 'top_y': 122.09, 'top_z': 375.95001},
 112: {'top_x': 53.019001, 'top_y': 122.052, 'top_z': 375.81},
 113: {'top_x': -23.993999, 'top_y': 109.937, 'top_z': 376.48999},
 114: {'top_x': -9.9799995, 'top_y': 109.944, 'top_z': 376.32999},
 115: {'top_x': 3.9960001, 'top_y': 109.951, 'top_z': 376.22},
 116: {'top_x': 18.025, 'top_y': 109.959, 'top_z': 376.07001},
 117: {'top_x': 32.023998, 'top_y': 109.944, 'top_z': 375.95999},
 118: {'top_x': 46.021999, 'top_y': 109.949, 'top_z': 375.81},
 119: {'top_x': -16.987, 'top_y': 97.809998, 'top_z': 376.35999},
 120: {'top_x': -2.9979999, 'top_y': 97.819, 'top_z': 376.23001},
 121: {'top_x': 11.019, 'top_y': 97.833, 'top_z': 376.14001},
 122: {'top_x': 25.004, 'top_y': 97.835999, 'top_z': 376.01001},
 123: {'top_x': 39.021, 'top_y': 97.82, 'top_z': 375.89999},
 124: {'top_x': -9.9870005, 'top_y': 85.702003, 'top_z': 376.25},
 125: {'top_x': 4.0180001, 'top_y': 85.697998, 'top_z': 376.16},
 126: {'top_x': 18.018, 'top_y': 85.716003, 'top_z': 376.04001},
 127: {'top_x': 32.013, 'top_y': 85.712997, 'top_z': 375.91}
}

EastHex = np.array([[ 57, 58, 59, 60, -1, -1, -1],
                    [ 61, 62, 63, 64, 65, -1, -1],
                    [ 66, 67, 68, 69, 70, 71, -1],
                    [ 72, 73, 74, -1, 75, 76, 77],
                    [ -1, 78, 79, 80, 81, 82, 83],
                    [ -1, -1, 84, 85, 86, 87, 88],
                    [ -1, -1, -1, 89, 90, 91, 92]])

SouthHex = np.array([[ -1, 93, 94, 95, -1, -1, -1],
                     [ 96, 97, 98, 99,100, -1, -1],
                     [101,102,103,104,105,106, -1],
                     [107,108,109, -1,110,111,112],
                     [ -1,113,114,115,116,117,118],
                     [ -1, -1,119,120,121,122,123],
                     [ -1, -1, -1,124,125,126,127]])

tile_info = {
 0: {'Tile': 11, 'cable': 90},
 1: {'Tile': 12, 'cable': 90},
 2: {'Tile': 13, 'cable': 150},
 3: {'Tile': 14, 'cable': 150},
 4: {'Tile': 15, 'cable': 150},
 5: {'Tile': 16, 'cable': 90},
 6: {'Tile': 17, 'cable': 150},
 7: {'Tile': 18, 'cable': 90},
 8: {'Tile': 21, 'cable': 230},
 9: {'Tile': 22, 'cable': 230},
 10: {'Tile': 23, 'cable': 230},
 11: {'Tile': 24, 'cable': 230},
 12: {'Tile': 25, 'cable': 230},
 13: {'Tile': 26, 'cable': 230},
 14: {'Tile': 27, 'cable': 150},
 15: {'Tile': 28, 'cable': 230},
 16: {'Tile': 31, 'cable': 230},
 17: {'Tile': 32, 'cable': 230},
 18: {'Tile': 33, 'cable': 150},
 19: {'Tile': 34, 'cable': 150},
 20: {'Tile': 35, 'cable': 150},
 21: {'Tile': 36, 'cable': 150},
 22: {'Tile': 37, 'cable': 150},
 23: {'Tile': 38, 'cable': 150},
 24: {'Tile': 41, 'cable': 150},
 25: {'Tile': 42, 'cable': 150},
 26: {'Tile': 43, 'cable': 150},
 27: {'Tile': 44, 'cable': 150},
 28: {'Tile': 45, 'cable': 150},
 29: {'Tile': 46, 'cable': 230},
 30: {'Tile': 47, 'cable': 320},
 31: {'Tile': 48, 'cable': 230},
 32: {'Tile': 61, 'cable': 150},
 33: {'Tile': 62, 'cable': 150},
 34: {'Tile': 63, 'cable': 150},
 35: {'Tile': 64, 'cable': 150},
 36: {'Tile': 65, 'cable': 230},
 37: {'Tile': 66, 'cable': 150},
 38: {'Tile': 67, 'cable': 150},
 39: {'Tile': 68, 'cable': 150},
 40: {'Tile': 81, 'cable': 230},
 41: {'Tile': 82, 'cable': 150},
 42: {'Tile': 83, 'cable': 230},
 43: {'Tile': 84, 'cable': 230},
 44: {'Tile': 85, 'cable': 230},
 45: {'Tile': 86, 'cable': 230},
 46: {'Tile': 87, 'cable': 230},
 47: {'Tile': 88, 'cable': 90},
 48: {'Tile': 91, 'cable': 320},
 49: {'Tile': 92, 'cable': 230},
 50: {'Tile': 93, 'cable': 150},
 51: {'Tile': 94, 'cable': 230},
 52: {'Tile': 95, 'cable': 230},
 53: {'Tile': 96, 'cable': 230},
 54: {'Tile': 97, 'cable': 230},
 55: {'Tile': 98, 'cable': 230},
 56: {'Tile': 99, 'cable': 90},
 57: {'Tile': 1001, 'cable': 90},
 58: {'Tile': 1002, 'cable': 90},
 59: {'Tile': 1003, 'cable': 90},
 60: {'Tile': 1004, 'cable': 150},
 61: {'Tile': 1005, 'cable': 150},
 62: {'Tile': 1006, 'cable': 150},
 63: {'Tile': 1007, 'cable': 150},
 64: {'Tile': 1008, 'cable': 150},
 65: {'Tile': 1009, 'cable': 150},
 66: {'Tile': 1010, 'cable': 90},
 67: {'Tile': 1011, 'cable': 150},
 68: {'Tile': 1012, 'cable': 150},
 69: {'Tile': 1013, 'cable': 150},
 70: {'Tile': 1014, 'cable': 150},
 71: {'Tile': 1015, 'cable': 150},
 72: {'Tile': 1016, 'cable': 90},
 73: {'Tile': 1017, 'cable': 90},
 74: {'Tile': 1018, 'cable': 150},
 75: {'Tile': 1019, 'cable': 150},
 76: {'Tile': 1020, 'cable': 150},
 77: {'Tile': 1021, 'cable': 230},
 78: {'Tile': 1022, 'cable': 90},
 79: {'Tile': 1023, 'cable': 150},
 80: {'Tile': 1024, 'cable': 150},
 81: {'Tile': 1025, 'cable': 150},
 82: {'Tile': 1026, 'cable': 150},
 83: {'Tile': 1027, 'cable': 150},
 84: {'Tile': 1028, 'cable': 90},
 85: {'Tile': 1029, 'cable': 150},
 86: {'Tile': 1030, 'cable': 150},
 87: {'Tile': 1031, 'cable': 150},
 88: {'Tile': 1032, 'cable': 150},
 89: {'Tile': 1033, 'cable': 150},
 90: {'Tile': 1034, 'cable': 150},
 91: {'Tile': 1035, 'cable': 150},
 92: {'Tile': 1036, 'cable': 150},
 93: {'Tile': 1038, 'cable': 90},
 94: {'Tile': 1039, 'cable': 90},
 95: {'Tile': 1040, 'cable': 150},
 96: {'Tile': 1041, 'cable': 150},
 97: {'Tile': 1042, 'cable': 150},
 98: {'Tile': 1043, 'cable': 150},
 99: {'Tile': 1044, 'cable': 150},
 100: {'Tile': 1045, 'cable': 150},
 101: {'Tile': 1046, 'cable': 90},
 102: {'Tile': 1047, 'cable': 150},
 103: {'Tile': 1048, 'cable': 150},
 104: {'Tile': 1049, 'cable': 150},
 105: {'Tile': 1050, 'cable': 150},
 106: {'Tile': 1051, 'cable': 230},
 107: {'Tile': 1052, 'cable': 90},
 108: {'Tile': 1053, 'cable': 90},
 109: {'Tile': 1054, 'cable': 150},
 110: {'Tile': 1055, 'cable': 150},
 111: {'Tile': 1056, 'cable': 150},
 112: {'Tile': 1057, 'cable': 230},
 113: {'Tile': 1058, 'cable': 90},
 114: {'Tile': 1059, 'cable': 150},
 115: {'Tile': 1060, 'cable': 150},
 116: {'Tile': 1061, 'cable': 150},
 117: {'Tile': 1062, 'cable': 150},
 118: {'Tile': 1063, 'cable': 230},
 119: {'Tile': 1064, 'cable': 90},
 120: {'Tile': 1065, 'cable': 150},
 121: {'Tile': 1066, 'cable': 150},
 122: {'Tile': 1067, 'cable': 150},
 123: {'Tile': 1068, 'cable': 150},
 124: {'Tile': 1069, 'cable': 150},
 125: {'Tile': 1070, 'cable': 150},
 126: {'Tile': 1071, 'cable': 150},
 127: {'Tile': 1072, 'cable': 150},
}
