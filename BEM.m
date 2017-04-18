function [thrust torque power] = BEM(liftfunc, chord, pitch, roottwist, tiptwist, BLADE, V, RPM, rho, blades)
    
    % build linear twist function and dr
    element = BLADE(end)-BLADE(end-1)
    omega = (RPM/60)*2*pi;
    m = (roottwist - tiptwist)/(BLADE(1)-BLADE(end));
    const = roottwist-m*BLADE(1);
    
    % calculate equation constants outside while-loop
    K1 = 0.5*rho*blades*chord;
    K2 = 0.5*rho*blades*chord*BLADE;
    K3 = 4.0*pi*rho*V.^2;
    K4 = K3*omega./V;   
    
    % allocating memory
    A=0.1*ones(size(BLADE));
    B=0.01*ones(size(BLADE));
    ANEW=zeros(size(BLADE));
    BNEW=zeros(size(BLADE));
    DtDr=zeros(size(BLADE));
    DqDr=zeros(size(BLADE));
    
    % initialize thrust and torque vars, relaxation factor, and control variables
    thrust = 0;
    torque = 0;
    convfact = 0.1;    
    finished=0;
    count = 0;
    
    % calculate theta along the blade, taking account of linear twist
    THETA = (m*BLADE+const+pitch)*pi/180;
    
    % calculate induced velocity and swirl
    while finished==0,
      V0=V*(1+A);
      V2=omega*BLADE.*(1-B);
      PHI=atan2(V0,V2);
      ALPHA=THETA-PHI;
      [CL CD] = liftfunc(ALPHA);
      VLOCALSQU=V0.^2+V2.^2;
      DtDr=K1*VLOCALSQU.*(CL.*cos(PHI)-CD.*sin(PHI));
      DqDr=K2.*VLOCALSQU.*(CD.*cos(PHI)+CL.*sin(PHI));
      TEM1=DtDr./(K3*BLADE.*(1+A));
      TEM2=DqDr./(K4*BLADE.^3.*(1+A));
      ANEW=convfact*A+(1-convfact)*TEM1;
      BNEW=convfact*B+(1-convfact)*TEM2;
      
      if (max(abs(ANEW-A))<1.0e-5),
        if (max(abs(BNEW-B))<1.0e-5),
          finished=1;
        end;
      end;
      
      A=ANEW;
      B=BNEW;
      count=count+1;
      
      if (count>500),
        finished=1;
      end;  
    end
    
    thrust = sum(DtDr*element);
    torque = sum(DqDr*element);
    power = torque*omega;
end