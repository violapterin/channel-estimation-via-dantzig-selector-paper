# Adjusted part

# # Set 1 (Uncomment to try these values)
# nn_yy = 3
# nn_rr = 6
# nn_hh = 9
# num_repeat = 12

# Set 2 (Uncomment to try these values)
nn_yy = 4
nn_rr = 8
nn_hh = 12
num_repeat = 16

# # Set 3 (Uncomment to try these values)
# nn_yy = 6
# nn_rr = 12
# nn_hh = 18
# num_repeat = 20

# # Set 4 (Uncomment to try these values)
# nn_yy = 8
# nn_rr = 16
# nn_hh = 24
# num_repeat = 24

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Fixed part

nn_y = nn_yy * nn_yy
nn_h = nn_hh * nn_hh
num_grid_phase = 16
lG_ant = 3
d_ant = 2
ll = 4
num_try_sigma = 2 * 3 + 1
num_try_gamma = 2 * 1 + 1 # DS
max_iter_oommpp = 4 * nn_h

