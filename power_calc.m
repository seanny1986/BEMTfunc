clear all;

rad = 0.085;
chord = 0.012;
stations = 10;
element = rad/stations;
BLADE = element/2:element:(rad-element/2);
BETA = 15*ones(size(BLADE));
CHORD = chord*ones(size(BLADE));
liftfunc = @NACA0012;
pitch = 0;
rpm_max = 15000;
rho = 1.225;
rho_isa = 1.225;
blades = 2;
max_discharge = 0.85;
avionics = 10;

dt = 0.1;
mass = 1.2;
gravity = 9.81;
motors = 4;
E = 6600;                                                                                     % mAh rating of battery
V = 14.8;                                                                                     % voltage rating of battery
hover_eta = 0.4;                                                                              % prop hover efficiency (can calc later)
max_alt = 500;
cd = 1.35;
Aeff = 0.11;
dia = 2*rad;
disk_area = pi*rad^2;                                                                         % disk area of a single rotor

no_runs = 5;

figure(1);
ax1 = gca;

figure(2);
ax2 = gca;

figure(3);
ax3 = gca;

figure(4);
ax4 = gca;

figure(5);
ax5 = gca;

C = {'k','b','r','g','m'};

for i=0:(no_runs-1)
  t = 0;
  v = 0.0001;
  x = 0.0001;
  TIME = [];
  XS = [];
  VS = [];
  THRUST = [];
  POWER = [];
  EFFICIENCY = [];
  rpm = rpm_max-i*250;
  n = rpm/60;
  while x<max_alt
    [thrust torque power] = BEM(liftfunc, pitch, CHORD, BETA, BLADE, v, rpm, rho, blades);
    thr = 4*thrust;
    DRAG = Aeff*cd*0.5*rho*v^2;
    dvdt = (thr-DRAG-mass*gravity)/mass;
    v = v+dvdt*dt;
    x = x+v*dt;
    XS = [XS, x];
    VS = [VS, v];
    pwr = 4*power+mass*gravity*v+DRAG*v+avionics;
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
  v
  4*power
  total_energy = sum(POWER)*dt
  hover_power = (motors/4/hover_eta)*rho*disk_area*(2*mass*gravity/motors/rho/disk_area)^(3/2)      % momentum theory hover power
  hover_time = (max_discharge*E*V*3.6-total_energy)/(hover_power+avionics)/60                                            % calc hover time in minutes
  
  figure(1);
  plot(ax1,TIME,XS,'color',C{i+1});
  hold on;
  
  figure(2);
  plot(ax2,TIME,VS,'color',C{i+1});
  hold on;

  figure(3);
  plot(ax3,TIME,THRUST,'color',C{i+1});
  hold on;
  
  figure(4);
  plot(ax4,TIME,POWER,'color',C{i+1});
  hold on;
  
  figure(5);
  plot(ax5,TIME,EFFICIENCY,'color',C{i+1});
  hold on;
end

figure(1);
title('Altitude vs Time');
xlabel('Time (s)');
ylabel('Altitude (m)');
xlim([0,t]);
legend('15000 rpm', '14750 rpm', '14500 rpm', '14250 rpm', '14000 rpm');
grid on;

figure(2);
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
xlim([0,t]);
legend('15000 rpm', '14750 rpm', '14500 rpm', '14250 rpm', '14000 rpm');
grid on;

figure(3);
title('Thrust vs Time');
xlabel('Time(s)');
ylabel('Thrust (N)');
xlim([0,t]);
ylim([0,50])
legend('15000 rpm', '14750 rpm', '14500 rpm', '14250 rpm', '14000 rpm');
grid on;

figure(4);
title('Power Draw vs Time');
xlabel('Time(s)');
ylabel('Power Draw (Watts)');
xlim([0,t]);
legend('15000 rpm', '14750 rpm', '14500 rpm', '14250 rpm', '14000 rpm');
grid on;

figure(5);
title('Efficiency vs Time');
xlabel('Time(s)');
ylabel('Efficiency (%)');
xlim([0,t]);
legend('15000 rpm', '14750 rpm', '14500 rpm', '14250 rpm', '14000 rpm');
grid on;