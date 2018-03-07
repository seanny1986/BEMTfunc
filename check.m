clear all;

rad = 0.085;
rho = 1.225;
rho_isa = 1.225;

dt = 0.1;
mass = 1.4;
gravity = 9.81;
motors = 4;
E = 6600;                                                                                     % mAh rating of battery
V = 14.8;                                                                                     % voltage rating of battery
hover_eta = 0.4;                                                                              % prop hover efficiency (can calc later)
max_alt = 1000;
cd = 1.35;
Aeff = 0.01;
dia = 2*rad;
n = rpm/60;
disk_area = pi*rad^2;                                                                         % disk area of a single rotor

no_runs = 5;

figure(1);
ax1 = gca;

t = 0;
v = 0.0001;
x = 0.0001;
TIME = [];
XS = [];
VS = [];
THRUST = [];
POWER = [];
EFFICIENCY = [];
rpm = 21189-i*1000;

T = 2.6/4;

while x<max_alt
  DRAG = Aeff*cd*0.5*rho*v^2;
  thr = T-DRAG-mass*gravity;
  dvdt = thr/mass;
  v = v+dvdt*dt;
  x = x+v*dt;
  XS = [XS, x];
  VS = [VS, v];
  pwr = climb_power(thr,v,rho,disk_area);
  THRUST = [THRUST, thr];
  POWER = [POWER, pwr];
  TIME = [TIME, t];
  J = v/n/dia;
  kt=thrust/(rho*n*n*dia*dia*dia*dia);
  kq=torque/(rho*n*n*dia*dia*dia*dia*dia);
  EFFICIENCY=[EFFICIENCY, (J/2/pi)*kt/kq];     
  t = t+dt;
  isa = atmosisa(x);
  rho = isa(4)*rho_isa;
end