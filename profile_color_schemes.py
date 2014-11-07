# ====================================
# ADDON TO CHANGE POFILEFACE COLOR:
# ====================================
import colorsys
from PyQt4 import QtGui

# Pick two colors. Values from 0 to 1. See "hue" at
# http://en.wikipedia.org/wiki/HSL_and_HSV
COLOR1= 0.1
COLOR2= 0.5
COLOR_INTENSITY = 0.6
def gradient(hue):
    min_lightness = 0.35 
    max_lightness = 0.85
    base_value = COLOR_INTENSITY

    # each gradient must contain 100 lightly descendant colors
    colors = []   
    rgb2hex = lambda rgb: '#%02x%02x%02x' % rgb
    l_factor = (max_lightness-min_lightness) / 100.
    l = min_lightness
    while l<=max_lightness:
        l += l_factor
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), 
                                 colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
    return colors

def get_color_gradient(self):
    colors = []
    for c in  gradient(COLOR1):
        color=QtGui.QColor(c)
        colors.append(color)
    colors.append(QtGui.QColor("white"))
    for c in  reversed(gradient(COLOR2)):
        color=QtGui.QColor(c)
        colors.append(color)
    return colors 



# ====================================
# JUST A TEST SCRIPT
# ====================================

from ete2 import ClusterTree, ProfileFace

# Let's replace the function that generates the color gradients in
# ProfileFaces, so the config is applied in all profile faces. 
ProfileFace.get_color_gradient = get_color_gradient

# Test it with a clustering tree!
matrix = """
#Names\tcol1\tcol2\tcol3\tcol4\tcol5\tcol6\tcol7
A\t-1.23\t-0.81\t1.79\t0.78\t-0.42\t-0.69\t0.58
B\t-1.76\t-0.94\t1.16\t0.36\t0.41\t-0.35\t1.12
C\t-2.19\t0.13\t0.65\t-0.51\t0.52\t1.04\t0.36
D\t-1.22\t-0.98\t0.79\t-0.76\t-0.29\t1.54\t0.93
E\t-1.47\t-0.83\t0.85\t0.07\t-0.81\t1.53\t0.65
F\t-1.04\t-1.11\t0.87\t-0.14\t-0.80\t1.74\t0.48
G\t-1.57\t-1.17\t1.29\t0.23\t-0.20\t1.17\t0.26
H\t-1.53\t-1.25\t0.59\t-0.30\t0.32\t1.41\t0.77
"""
t = ClusterTree("(((A,B),(C,(D,E))),(F,(G,H)));", text_array=matrix)
t.show("heatmap")
