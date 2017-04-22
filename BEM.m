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
    ERRA=ones(size(BLADE)); ERRB=ones(size(BLADE));
    
    % specify step size for gradient descent, error threshold
    stepa = 0.001*ones(size(BLADE)); stepb=0.0001*ones(size(BLADE)); thresh = 1e-7;
    
    % use gradient descent to find the solution
    function [A0 B0 ERRA ERRB] = bemfuncmin(A0, B0, stepa, stepb)
      V0=V*(1+A0);  
      V2=omega*BLADE.*(1-B0);
      PHI=atan2(V0,V2);
      ALPHA=THETA-PHI;
      [CL CD]=liftfunc(ALPHA);
      VLOCSQU=V0.^2+V2.^2;
      DtDr=K1*VLOCSQU.*(CL.*cos(PHI)-CD.*sin(PHI));
      DqDr=K2.*VLOCSQU.*(CD.*cos(PHI)+CL.*sin(PHI));
      TEM1=DtDr./(K3*BLADE.*(1+A0));
      TEM2=DqDr./(K4*BLADE.^3.*(1+A0));
      dERRA_dA0=2*(A0-TEM1).*(1+DtDr./K3./BLADE./(1+A0).^2);
      dERRB_dB0=2*(B0-TEM2);
      ERRA=(A0-TEM1).^2;
      ERRB=(B0-TEM2).^2;
      A0=A0-dERRA_dA0.*stepa;
      B0=B0-dERRB_dB0.*stepb;
    end
    
    % solve for a known A and B value
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
    
    % gradient descent
    while sum(ERRA.^2)>=thresh^2 & sum(ERRB.^2)>=thresh^2,
      [A0 B0 ERRA ERRB]=bemfuncmin(A0, B0, stepa, stepb);
    end
    
    [DT DQ] = bemsolve(A0, B0);
    thrust=sum(DT*element);
    torque=sum(DQ*element);
    power = torque*omega;
end
