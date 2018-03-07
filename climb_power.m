function p = climb_power(thrust, velocity, rho, area)
  p = thrust*(velocity/2+sqrt((velocity/2)^2+thrust/2/rho/area))
end