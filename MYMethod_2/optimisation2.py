import numpy as np
def objective(parameters, x_coords, y_coords, x_coords_pic, y_coords_pic):
    
    mod_principal_x, mod_principal_y, mod_focal_len_x, mod_focal_len_y, rad_dist_k2, rad_dist_k4, rad_dist_k6, tang_dist_p1, tang_dist_p2 = parameters

    x_mod = (mod_focal_len_x * x_coords_pic) + mod_principal_x 
    y_mod = (mod_focal_len_y * x_coords_pic) + mod_principal_y

    r = pow((pow(x_mod, 2) + pow(y_mod, 2)), 0.5)

    x_cal = (x_mod * (1 + rad_dist_k2 * (pow(r, 2)) + rad_dist_k4 * (pow(r, 2)) + rad_dist_k6 * (pow(r, 2)))) + (2 * tang_dist_p1 * x_mod * y_mod + tang_dist_p2 * (pow(r, 2) + 2 * pow(x_mod,2)))
    y_cal = (x_mod * (1 + rad_dist_k2 * (pow(r, 2)) + rad_dist_k4 * (pow(r, 2)) + rad_dist_k6 * (pow(r, 2)))) + (2 * tang_dist_p2 * x_mod * y_mod + tang_dist_p1 * (pow(r, 2) + 2 * pow(y_mod,2)))

    print(x_cal, y_cal)

    errors = pow((pow((x_coords - x_cal), 2) + pow((y_coords - y_cal), 2)), 0.5)

    return errors