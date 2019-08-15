import numpy as np

# # Part 1

# precoder dimension
nn_yy = 3
nn_rr = 6
nn_hh = 9
nn_y = nn_yy *nn_yy
nn_h = nn_hh *nn_hh
num_grid_phase = 16

# channel parameters
lG_ant = 3
d_ant = 2
ll = 4

# simulation count
num_scale_noise = 9
num_try_gG = 3
num_repeat = 12

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Part 2
# # Uncomment the following to try different values:

# # precoder dimension
# nn_yy = 4
# nn_rr = 8
# nn_hh = 12
# nn_y = nn_yy *nn_yy
# nn_h = nn_hh *nn_hh
# num_grid_phase = 16

# # channel parameters
# lG_ant = 3
# d_ant = 2
# ll = 4

# # simulation count
# num_scale_noise = 12
# num_try_gG = 4
# num_repeat = 18

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Part 3
# # Uncomment the following to try different values:

# # precoder dimension
# nn_yy = 8
# nn_rr = 12
# nn_hh = 16
# nn_y = nn_yy *nn_yy
# nn_h = nn_hh *nn_hh
# num_grid_phase = 16

# # channel parameters
# lG_ant = 3
# d_ant = 2
# ll = 4

# # simulation count
# num_scale_noise = 15
# num_try_gG = 6
# num_repeat = 24
