function [thrust torque power] = BEM(liftfunc, pitch, CHORD, BETA, BLADE, v, rpm, rho, blades)
    % args -- liftfunc is a function handle to a lookup table that takes in an angle, and produces [CL, CD]
    % pitch is the pitching angle of the blade in degrees
    % CHORD is the vector of station-wise chords along the blade
    % BETA is the vector of station-wise blade angles due to twist, in degrees
    % BLADE is the vector of radial positions of each blade element
    % v is the entry velocity to the disc, in m/s
    % rpm is the rotation per minute of the prop_test_2
    % rho is atmospheric density in metric units
    % blades is the number of prop blades
    
    element = BLADE(end)-BLADE(end-1); omega = (rpm/60)*2*pi;                   % calculate element size and angular velocity
    K1 = 0.5*rho*blades*CHORD; K2 = 0.5*rho*blades*CHORD.*BLADE;                % calculate equation constants
    K3 = 4.0*pi*rho*v^2; K4 = K3*omega/v;                                       % calculate equation constants
    THETA = (pitch+BETA)*pi/180;                                                % calculate theta along the blade
    A0 = 0.1*ones(size(BLADE)); B0 = 0.01*ones(size(BLADE));                    % initialize A0, B0
    
    % calculate cost function and gradient at each element along the blade
    function [J grad] = bemfuncmin(ABs)
      A0 = ABs(1:size(BLADE,2)); B0 = ABs((size(BLADE,2)+1):end);               % unpack the A0 and B0 vectors
      A0 = A0'; B0 = B0';                                                       % put into row vector form
      [DtDr DqDr]=bemsolve(A0, B0);                                             % solve system of equations using estimated A and B values
      TEM1=DtDr./(K3*BLADE.*(1+A0)); TEM2=DqDr./(K4*BLADE.^3.*(1+A0));          % calculated A and B values from system of equations
      ERRA = 0.5*(A0-TEM1).^2; ERRB = 0.5*(B0-TEM2).^2;                         % error function
      dERRA_dA0 = (A0-TEM1).*(1+DtDr./K3./BLADE./(1+A0).^2);                    % derivative of the squared error function
      dERRB_dB0 = (B0-TEM2);                                                    % derivative of the squared error function
      J = sum([ERRA(:); ERRB(:)]); grad = [dERRA_dA0(:); dERRB_dB0(:)];         % convert back to column vectors
    end
    
    % solve lift and torque vectors for given A and B vectors
    function [DtDr DqDr] = bemsolve(A, B)
      V0=v*(1+A); V2=omega*BLADE.*(1-B); VLOCSQU=V0.^2+V2.^2;                   % calculate velocities
      PHI=atan2(V0,V2); ALPHA=THETA-PHI;                                        % calculate angles
      [CL CD]=liftfunc(ALPHA);                                                  % get lift and drag coefficients from our lift function
      DtDr=K1.*VLOCSQU.*(CL.*cos(PHI)-CD.*sin(PHI));                            % calculate thrust per element
      DqDr=K2.*VLOCSQU.*(CD.*cos(PHI)+CL.*sin(PHI));                            % calculate torque per element
    end
    
    costfunction = @(P)bemfuncmin(P);                                           % establish cost function
    options = optimset('MaxIter', 100000, 'GradObj','on');                      % set options
    ABs = [A0(:); B0(:)];                                                       % pack A0, B0 into single column vector
    [ABs cost] = fminunc(costfunction, ABs, options);                           % solving using unconstrained fmin     
    A0 = ABs(1:size(BLADE,2)); B0 = ABs((size(BLADE,2)+1):end);                 % unpack the A0 and B0 vectors
    A0 = A0'; B0 = B0';                                                         % put into row vector form
    [DT DQ] = bemsolve(A0, B0);                                                 % solve system using converged A0 and B0 values
    thrust=sum(DT*element); torque=sum(DQ*element);                             % calculate thrust, torque 
    power = torque*omega;                                                       % calculate power
end