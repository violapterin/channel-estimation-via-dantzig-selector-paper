import matplotlib.pyplot as plt # plotting functions
import numpy as np
import os # getcwd

# function draw (lst_x, lst_lst_y, str_label_x, lst_str_label_y, str_title)
len_line =size(lst_lst_y)

# to set axis & title
plt.close("all")
fig =plt.figure()
#plt.title( str_title, fontsize =15 )
#plt.xlabel( str_label_x, fontsize =12 )
#plt.ylabel( lst_str_label_y, fontsize =12 )

# to set plot style
num_style =4
# from the thinner to thicker
lst_style =['-', '--', '-.', ':']
num_color =5
# red, green, cyan, blue, black
lst_color =['r', 'g', 'c', 'b', 'k']
num_marker =6
# triangle down, triangle up, circle, square, star, diamond
lst_marker =['v', '^', 'o', 's', '*', 'D']
size_marker =7 # fixed marker size
width_line =3 # fixed line width

arr_x =np.array (arr_x)

for i in range (len_line):
    arr_y_i =np.array(lst_lst_y(i))
    plt.plot (arr_x,
              arr_y_i,
              markersize =size_marker,
              linewidth =width_line,
              linestyle =list_line_style[ int( idx % num_line_style ) ],
              color =list_color[ int( idx % num_color ) ],
              marker =list_marker[ int( idx % num_marker ) ],
              )

out_fig_path =os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))
              +"/plt/" +str_title +".png"
#print (out_fig_path)
plt.savefig (out_fig_path, bbox_inches ="tight")


