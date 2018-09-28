#!/usr/bin/env python
# the shebang line indicates this file be passed as a python3 script

import matplotlib.pyplot as plt # plotting functions
import os # getcwd

# to set input & output paths
top_dir = os.getcwd() # current working directory
out_data_path =top_dir +"/dat/table.txt"
out_fig_path =top_dir +"/plt/fig.png"
#print(out_fig_path)

# to retrieve line by line
handle_text =open(out_data_path)
content =handle_text.read()
list_line =content.split('\n') # newline as delimiter

# to set axis & title
#plt.close("all")
fig =plt.figure()
plt.title( "SA", fontsize =15 )
plt.xlabel( "Time", fontsize =12 )
plt.ylabel( "Rate", fontsize =12 )

# to set plot style
num_line_style =4
list_line_style =['-', '--', '-.', ':'] # the larger, the thicker
num_color =5
list_color =['r', 'g', 'c', 'b', 'k'] # red, green, cyan, blue, black
num_marker =6
list_marker =['v', '^', 'o', 's', '*', 'D'] # triangle down, triangle up, circle, square, star, diamond
size_width =3 # fixed line width
size_marker =7 # fixed marker size

vec_x =[] # to set empty list to enlarge its scope
# split w/o arguments parses any number of spaces and tabs
list_str =list_line[0]. split()
for tmp_str in list_str:
    # strip decimal mark to determine, whether a decimal number
    if not tmp_str. replace('-',''). replace('.',''). isdigit():
        continue
    else:
        vec_x.append(float(tmp_str))
# end for
list_line.pop(0)

for idx, line in enumerate(list_line):
    list_str =line.split()
    if not list_str: # "list_str" is empty
        continue
    # end if empty

    vec_y =[] # to create list to enlarge its scope
    
    for tmp_str in list_str:
        # to strip decimal mark & minus sign to check whether a decimal number
        if not (tmp_str. replace('-',''). replace('.',''). isdigit()):
            continue
        else:
            vec_y.append(float(tmp_str))
    # end for each word

    plt.plot( vec_x, vec_y,
            markersize =size_marker,
            linewidth =size_width,
            linestyle =list_line_style[ int( idx % num_line_style ) ],
            color =list_color[ int( idx % num_color ) ],
            marker =list_marker[ int( idx % num_marker ) ],
            )
# end for each line

# to set legend and save figure
#distShift =-9.0 # shifted, to the left, considerably
#plt.legend( loc='center left', borderaxespad =distShift )
#fig.show()
fig.savefig( out_fig_path, bbox_inches ="tight" )
#plt.show()

