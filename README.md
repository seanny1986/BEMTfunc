# BEMTfunc
Vectorized BEMT function for MATLAB/Octave. 

Takes a function handle for the lift and drag coefficient data; this should be of the form [cl cd]=coeff(alpha), and is passed to the BEMT function as @coeff. The BEMT function will call @coeff using an angle in radians. 

The BLADE input variable should be a vector of the form rstart:rstep:rstop. These should correspond to the midpoints of each element. For example, an 80cm blade with 10 elements, starting from 0 and ending at 0.8m, would be the vector 0.04:0.08:0.76.

The BEMT function will iterate to find the induced velocity at each point along the blade (this is computed for all elements simultaneously). Gradient descent is used to converge on a solution, and this solution is then used to calculate the thrust, torque, and power of the propeller.
