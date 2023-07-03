import numpy as np
def objective(parameters, x_coords, y_coords, x_ideal, y_ideal, focal_length):
    
    off_x, off_y, mod_fl, k_2, k_4, k_6, delta_phi, delta_omega, delta_kappa = parameters

    x_mod = x_coords + off_x
    y_mod = y_coords + off_y

    r = pow((pow(x_mod, 2) + pow(y_mod, 2)),0.5)
    x_ideal_cal = x_mod * (1 + (k_2 * pow(r, 2)) + (k_4 * pow(r, 4)) + (k_6 * pow(r, 6)))
    y_ideal_cal = y_mod * (1 + (k_2 * pow(r, 2)) + (k_4 * pow(r, 4)) + (k_6 * pow(r, 6)))

    rotation_matrix = np.array([[np.cos(delta_kappa)*np.cos(delta_phi) - np.sin(delta_kappa)*np.sin(delta_omega) + np.cos(delta_kappa)*np.sin(delta_phi)*np.sin(delta_omega), 
                             np.sin(delta_kappa)*np.cos(delta_phi) + np.cos(delta_kappa)*np.cos(delta_omega)*np.sin(delta_phi) - np.sin(delta_kappa)*np.sin(delta_phi)*np.sin(delta_omega),
                             -np.sin(delta_phi)*np.cos(delta_phi)*np.sin(delta_omega) + np.cos(delta_phi)*np.cos(delta_omega)],
                             [-np.sin(delta_kappa)*np.cos(delta_phi) - np.cos(delta_kappa)*np.cos(delta_omega)*np.sin(delta_phi) + np.sin(delta_kappa)*np.sin(delta_phi)*np.sin(delta_omega),
                              np.cos(delta_kappa)*np.cos(delta_phi) + np.sin(delta_kappa)*np.cos(delta_omega) + np.cos(delta_kappa)*np.sin(delta_phi)*np.sin(delta_omega), np.sin(delta_phi)*np.cos(delta_phi)*np.cos(delta_omega) + np.sin(delta_omega)*np.cos(delta_phi)],
                              [np.sin(delta_phi)*np.cos(delta_omega), -np.sin(delta_phi)*np.sin(delta_omega), np.cos(delta_phi)]])
    initial_vector = np.array([x_ideal, y_ideal, focal_length])
    rotated_vector = rotation_matrix @ initial_vector
    x_ideal_rot = rotated_vector[0]
    y_ideal_rot = rotated_vector[1]
    focal_length_rot = rotated_vector[2]

    x_ideal_mod = x_ideal_rot * (mod_fl/focal_length_rot)
    y_ideal_mod = y_ideal_rot * (mod_fl/focal_length_rot)

    errors = pow(pow((x_ideal_cal - x_ideal_mod),2) + pow((y_ideal_cal - y_ideal_mod),2),0.5)

    rmse = np.sqrt(np.mean(pow(errors,2)))

    return rmse