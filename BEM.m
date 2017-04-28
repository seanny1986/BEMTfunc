function [thrust torque power] = BEM(liftfunc, chord, pitch, roottwist, tiptwist, BLADE, V, RPM, rho, blades)
    
    % build linear twist function and dr
    element = BLADE(end)-BLADE(end-1);
    omega = (RPM/60)*2*pi;
    m = (roottwist - tiptwist)/(BLADE(1)-BLADE(end));
    const = roottwist-m*BLADE(1);
    
    % calculate equation constants outside while-loop
    K1 = 0.5*rho*blades*chord;
    K2 = 0.5*rho*blades*chord*BLADE;
    K3 = 4.0*pi*rho*V.^2;
    K4 = K3*omega./V;   
    
    % calculate theta along the blade, taking account of linear twist
    THETA = (m*BLADE+const+pitch)*pi/180;
        
    % allocating memory
    A0 = 0.1*ones(size(BLADE)); B0 = 0.01*ones(size(BLADE));
    TEM1 = zeros(size(BLADE)); TEM2 = zeros(size(BLADE));
    ERRA=ones(size(BLADE)); ERRB=ones(size(BLADE));
    runA=ones(size(BLADE)); runB=ones(size(BLADE));
    
    % specify step size for gradient descent, specify error threshold
    stepa = 0.01*ones(size(BLADE)); stepb=0.01*ones(size(BLADE)); thresh = 1e-7;
    
    % use gradient descent to find the induced velocity at each element
    function [A0 B0 TEM1 TEM2] = bemfuncmin(A0, B0, stepa, stepb)
      
      % solve system of equations using estimated A and B values
      [DtDr DqDr]=bemsolve(A0, B0);
      
      % calculated A and B values from system of equations
      TEM1=DtDr./(K3*BLADE.*(1+A0));
      TEM2=DqDr./(K4*BLADE.^3.*(1+A0));
      
      % derivative of the squared error function
      dERRA_dA0=2*(A0-TEM1).*(1+DtDr./K3./BLADE./(1+A0).^2);
      dERRB_dB0=2*(B0-TEM2);
      
      % update values of A0 and B0
      A0=A0-dERRA_dA0.*stepa;
      B0=B0-dERRB_dB0.*stepb;
      
      % might want to take measures to adjust step size -- allows for faster
      % convergence, and minimizes risk of blow up if the chosen step size is
      % too large
    end
    
    % solve lift and torque vectors for given A and B vectors
    function [DtDr DqDr] = bemsolve(A, B)
      V0=V*(1+A);  
      V2=omega*BLADE.*(1-B);
      PHI=atan2(V0,V2);
      ALPHA=THETA-PHI;
      [CL CD]=liftfunc(ALPHA);
      VLOCSQU=V0.^2+V2.^2;
      DtDr=K1*VLOCSQU.*(CL.*cos(PHI)-CD.*sin(PHI));
      DqDr=K2.*VLOCSQU.*(CD.*cos(PHI)+CL.*sin(PHI));
    end
    
    % gradient descent using bemfuncmin
    while sum(runA)>0 & sum(runB)>0,
      [A0 B0 TEM1 TEM2]=bemfuncmin(A0, B0, stepa, stepb);
      runA = abs(A0-TEM1)>thresh;
      runB = abs(B0-TEM2)>thresh;
    end
    
    % use converged A and B values to calculate thrust, torque, power
    [DT DQ] = bemsolve(A0, B0);
    thrust=sum(DT*element);
    torque=sum(DQ*element);
    power = torque*omega;
end
