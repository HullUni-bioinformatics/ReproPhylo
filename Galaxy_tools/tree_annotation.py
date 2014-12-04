
#######################################################################
#                        COMMAND LINE PARSING
#######################################################################
import argparse as agp

desc = ("Annotate all the trees in the provided Project")
parser = agp.ArgumentParser(description=desc)

#----------------------------------------------------------------------
parser.add_argument("projfile",
                    help="A pickle file to read and write",
                    type=str,
                    default=None)
parser.add_argument("figs",
                    help="A directory to keep fig files",
                    type=str, 
                    default='.')
parser.add_argument("-a", "--altprojfile",
                    help="An alternative pickle file to write",
                    type=str, 
                    default=None)
parser.add_argument("-_", "--__",
                    help=("Turn off warnings for Galaxy"),
                    action="store_true")
parser.add_argument("--noninteractive",
                    help=("Overwrite reports without asking"),
                    action="store_true")
#----------------------------------------------------------------------
parser.add_argument("-m", "--rootmeta",
                    help="Metadata by which to dermine the outgroup",
                    type=str, 
                    default='mid')
parser.add_argument("-r", "--rootvalue",
                    help=("Within the metadata column, determine a value"+
                          " that specifies the outgroup"),
                    type=str, 
                    default='mid')
parser.add_argument("-l", "--leaflabelmeta",
                    help=("Metadata to use as labels"),
                    type=str, nargs='+',
                    default=['source_organism'])
parser.add_argument("-t", "--leafcolormeta",
                    help=("Metadata to use to determine label colours"),
                    type=str, 
                    default=None)
parser.add_argument("-s", "--leafcolors",
                    help=("contrast, blues, greens, reds, purples, browns"),
                    type=str, choices=['contrast',
                                       'blues',
                                       'greens',
                                       'reds',
                                       'purples',
                                       'browns'], 
                    default=None)
parser.add_argument("-g", "--nodebgmeta",
                    help=("Metadata to determine background colors"),
                    type=str,
                    default=None)
parser.add_argument("-i", "--nodebgcolors",
                    help=("contrast, blues, greens, reds, purples, browns"),
                    type=str, choices=['contrast',
                                       'blues',
                                       'greens',
                                       'reds',
                                       'purples',
                                       'browns'], 
                    default=None)
parser.add_argument("-n", "--nodesupport",
                    help=("Triplets of values: eg, black 100 90 gray 90 80"),
                    type=str, nargs='+',
                    default=['black', '100', '99',
                             'dimgray', '99', '75',
                             'silver', '75', '50'])
parser.add_argument("-b", "--heatmapmeta",
                    help=("A list of metadata column names"),
                    type=str, nargs='+',
                    default=None)
parser.add_argument("-d", "--heatmapcolor",
                    help=("Heat map color scheme (0,1 or 2)"),
                    type=str, choices=['0','1','2'], 
                    default="2")
parser.add_argument("-u", "--multifurc",
                    help=("Support value under which to collapse node"),
                    type=str, 
                    default=None)
parser.add_argument("-e", "--treescale",
                    help=("Tree width"),
                    type=str, 
                    default="1000")
#----------------------------------------------------------------------
report = parser.add_mutually_exclusive_group()
report.add_argument("-v", "--verbose_report",
                   help=("Write full html report to filename"),
                   type=str,
                   default=None)
report.add_argument("--html",
                   help=("Write html of fig file links"),
                   type=str,
                   default=None)
#----------------------------------------------------------------------
parser.add_argument("-c", "--version_ctrl",
                    help=("perform version control in CWD"),
                    action="store_true")
#----------------------------------------------------------------------
cl = parser.parse_args()
#----------------------------------------------------------------------
# Coloring labels requires both metadata column and colors
if len([x for x in (cl.leafcolormeta,
                    cl.leafcolors) if x is not None])== 1:
    parser.error('To color labels, both --leafcolormeta and'+
                 ' --leafcolors are required')
#----------------------------------------------------------------------
# Coloring backgrounds requires both metadata column and colors
if len([x for x in (cl.nodebgmeta,
                    cl.nodebgcolors) if x is not None])== 1:
    parser.error('To color backgrounds, both --nodebgmeta and'+
                 ' --nodebgcolors are required')
#----------------------------------------------------------------------
if not len(cl.nodesupport)%3==0:
    parser.error('Node support annotation requires triplets consisting of: '+
                 'color top_range bottom_range')
#----------------------------------------------------------------------
from reprophylo import *

colors = {'greens': ['DarkOliveGreen',
                     'Olive',
                     'OliveDrab',
                     'YellowGreen',
                     'LimeGreen',
                     'Lime',
                     'LawnGreen',
                     'Chartreuse',
                     'GreenYellow',
                     'SpringGreen',
                     'MediumSpringGreen',
                     'LightGreen',
                     'PaleGreen',
                     'DarkSeaGreen',
                     'MediumSeaGreen',
                     'SeaGreen',
                     'ForestGreen',
                     'Green',
                     'DarkGreen'],
            'reds': ['LightSalmon',
                     'Salmon',
                     'DarkSalmon',
                     'LightCoral',
                     'IndianRed',
                     'Crimson',
                     'FireBrick',
                     'DarkRed',
                     'Red',
                     'OrangeRed',
                     'Tomato',
                     'Coral',
                     'DarkOrange',
                     'Orange'],
            'blues': ['LightSteelBlue',
                     'PowderBlue',
                     'LightBlue',
                     'SkyBlue',
                     'LightSkyBlue',
                     'DeepSkyBlue',
                     'DodgerBlue',
                     'CornflowerBlue',
                     'SteelBlue',
                     'RoyalBlue',
                     'Blue',
                     'MediumBlue',
                     'DarkBlue',
                     'Navy',
                     'MidnightBlue',
                     'MediumAquamarine',
                     'Aqua',
                     'Cyan',
                     'LightCyan',
                     'PaleTurquoise',
                     'Aquamarine',
                     'Turquoise',
                     'MediumTurquoise',
                     'DarkTurquoise',
                     'LightSeaGreen',
                     'CadetBlue',
                     'DarkCyan',
                     'Teal'],
            'purples': ['Lavender',
                     'Thistle',
                     'Plum',
                     'Violet',
                     'Orchid',
                     'Fuchsia',
                     'Magenta',
                     'MediumOrchid',
                     'MediumPurple',
                     'BlueViolet',
                     'DarkViolet',
                     'DarkOrchid',
                     'DarkMagenta',
                     'Purple',
                     'Indigo',
                     'DarkSlateBlue',
                     'SlateBlue',
                     'MediumSlateBlue'],
            'browns': ['Cornsilk',
                     'BlanchedAlmond',
                     'Bisque',
                     'NavajoWhite',
                     'Wheat',
                     'BurlyWood',
                     'Tan',
                     'RosyBrown',
                     'SandyBrown',
                     'Goldenrod',
                     'DarkGoldenrod',
                     'Peru',
                     'Chocolate',
                     'SaddleBrown',
                     'Sienna',
                     'Brown',
                     'Maroon'],
            'contrast': ['Red',
                         'Orange',
                         'Yellow',
                         'SaddleBrown',
                         'LimeGreen',
                         'Aqua',
                         'Blue',
                         'Fuchsia',
                         'Silver']
            
          }
    
def chekpoint_project(pj):
    if cl.altprojfile:
        pickle_pj(pj, cl.altprojfile)
    else:
        pickle_pj(pj, cl.projfile)

def colors_dict(pj, meta, colors_choice):
    dct = {}
    indx = 0
    for record in pj.records:
        for f in record.features:
            if meta in f.qualifiers.keys():
                if not f.qualifiers[meta][0] in dct.keys():
                    if indx < len(colors[colors_choice]):
                        dct[f.qualifiers[meta][0]] = colors[colors_choice][indx]
                        indx += 1
                    else:
                        dct[f.qualifiers[meta][0]] = colors[colors_choice][0]
                        indx = 1
    return dct
    

def dict_from_list_pairs(t):
    if t:
        dictionary = {}
        pairs = zip(t[::2], t[1::2])
        for i in pairs:
            dictionary[i[0]] = i[1]
        return dictionary
    else:
        return None
    
def dict_from_list_triplets(t):
    if t:
        dictionary = {}
        trios = zip(t[::3], t[1::3], t[2::3])
        for i in trios:
            dictionary[i[0]] = [int(i[1]),int(i[2])]
        return dictionary
    else:
        return None
    
def int_or_none(val):
    try:
        return int(val)
    except:
        return None

def removetype(path):
    import shutil
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)

#######################################################################
#                       Execution
#######################################################################
cline = ''

for i in sys.argv:
    cline += i+' '
print 'ReproPhylo was called with:'
print cline
print

if cl.__:
    warnings.simplefilter('ignore')

if cl.version_ctrl:
    start_git()
    
pj = unpickle_pj(cl.projfile)
pj.clear_tree_annotations()
pj.annotate(cl.figs,
            cl.rootmeta,
            cl.rootvalue,
            cl.leaflabelmeta,
            leaf_node_color_meta=cl.leafcolormeta,
            leaf_label_colors= colors_dict(pj, cl.leafcolormeta, cl.leafcolors),
            node_bg_meta=cl.nodebgmeta,
            node_bg_color=colors_dict(pj, cl.nodebgmeta, cl.nodebgcolors),
            node_support_dict=dict_from_list_triplets(cl.nodesupport),
            heat_map_meta=cl.heatmapmeta,
            heat_map_colour_scheme=cl.heatmapcolor,
            multifurc=int_or_none(cl.multifurc),
            scale=int(cl.treescale),
            html=cl.html
            )

if cl.verbose_report:
    if cl.noninteractive:
        if os.path.exists(cl.verbose_report):
            removetype(cl.verbose_report)
        if os.path.exists(cl.verbose_report.rpartition('.')[0]):
            removetype(cl.verbose_report.rpartition('.')[0])
        if os.path.exists(cl.verbose_report+'.zip'):
            removetype(cl.verbose_report+'.zip')
    publish(pj, cl.verbose_report, cl.figs)


chekpoint_project(pj)