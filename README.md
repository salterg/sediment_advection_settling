# sediment_advection_settling
Model solves the advection/settling equations in 2D for sediment with a van Leer flux limiter for advection.

To run the model, input hydrodynamic data is needed from the ANUGA model described by Wright et al. (2022).
The data should be saved as a .mat file, with:
x_full and y_full as row vectors of x and y coordinates (UTM in meters), respectively
t_full as column vector of model output times in seconds 
elev_full, h_full, u_full, v_full are ANUGA model outputs of elevation (m), depth (m), x velocity (m/s), y velocity (m/s)
structure of variables is each row represents a time from t_full, and each column corresponds to one of x,y points from x_full&y_full
