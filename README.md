# BEMTfunc
Blade Element Momentum Theory Function for MATLAB

Vectorized BEMT function for MATLAB/Octave. 

Takes a function handle for the lift and drag coefficient data; this should be of the form [cl cd]=coeff(alpha), and is passed to the BEMT function as @coeff. The BEMT function will call @coeff using an angle in radians. 

The BLADE input variable should be a vector of the form rstart:rstep:rstop. These correspond to the points should correspond to the midpoint of your element. For example, a 80cm blade with 10 elements would be the vector 0.04:0.08:0.76.

The BEMT function will iterate to find the induced velocity at each point along the blade. A relaxation value is used to help convergence; this is initially set to 0.1, but can easily be changed (look for the convfact variable). As is, it will iterate for a maximum of 500 iterations, though this can easily be changed in the code. There is no guarantee of convergence at any given element -- generally, smaller relaxation values with a higher number of maximum iterations will be more stable.
