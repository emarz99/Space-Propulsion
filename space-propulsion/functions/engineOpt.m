function [costVal] = engineOpt(x, engine)
    
    engine.OF = x(1);
    Tobj = engine.Thrust; 

    %%% Compute coefficients-----------------------------------------------
    engine.input.OF = x(1); 
    engine.input.subar = linspace(engine.CR, 1, engine.nSub+2);
    engine.input.subar = engine.input.subar(2:end-1); 
    
    coeff = createDataCEA(engine.input, engine.ceaPath); 
    p = coeff(1, :);    
    Temp = coeff(2, :); 
    rho = coeff(3, :);   
    gam = coeff(4, :); 
    vSon = coeff(5, :);
    Mach = coeff(6, :); 
    %kCond = coeff(7, :);
    %Pr = coeff(8, :);
    cStar = coeff(9, :); 
    cT = coeff(10, :);
    %Ivac = coeff(11, :);
    %Isp = coeff(12, :);
    %visc = coeff(13, :); 
     
    
    %% Computation -------------------------------------------------------
    
    % function to compute the throat area to obtain the given thrust
    AtFun = @(x) nozzleDesign(x, engine)*((engine.pc*x)/cStar(end))*Mach(end)*vSon(end) + p(end)*x*engine.eps - Tobj; 
    At0 = Tobj/(cT(end)*engine.pc);    % initial guesss throat area [m^2]
       
    [At, res] = fzero(AtFun, At0); 

    mDotP = (engine.pc*At)/cStar(end);         % propellant mass flow rate [kg/s]
    mDotF = mDotP/(engine.OF + 1);             % fuel mass flow rate [kg/s]
    mDotO = mDotF * engine.OF; 
    Ivac = Tobj/mDotP;                         % specific impulse

    %%% chamber design
    Ac = engine.CR * At;                        % chamber area [m2]
    Rc = sqrt(Ac/pi);                           % chamber radius [m]
    Vcc = engine.Lstar * At;                    % chamber + convergent volume [m3]

    l = Rc/tan(engine.thetac); 
    xC = Rc - sqrt(At/pi);                      % convergent length [m]
    Vconv = 1/3 * (Ac*l -  At*(l-xC));          % convergent volume [m3]

    Lc = (Vcc - Vconv)/Ac;                      % chamber length                       [m2]
    Aside = 2*Rc*pi*Lc;                         % compute the side area of the chamber [m2]

    % compute flow velocity in the chamber to satisfy the continuity
    % equation
    f = @(x) (mDotP/(Ac*rho(1)))*(1 + 0.5*(gam(1)-1)*(x/vSon(1))^2)^(1/(gam(1)-1)) - x; 
    
    vC = fzero(f, 0);                           % velocity inside the chamber [m/s]
    machC = vC/vSon(1);                         % mach inside the chamber     

    % obtain static temperature in the chamber
    [~, tSC, ~] = totalToStatic(gam(1), machC, p(1), Temp(1), rho(1));
    
    hC = hConvection(2*Rc, vC, coeff(:, 1), true);    % compute convection heat flux coefficient [W/(m2 K)]
    
    % Heat flux out from the combustion chamber [W]
    
    
    
    if isnan(engine.kCoat)
        QE = hC*Aside*(tSC(1) - engine.Twall); 
    else
        xIs = ((tSC - engine.Twall)/(tSC - engine.Tcoat) - 1) * engine.kCoat/hC;
    
        HC = 1/(1/hC + xIs/engine.kCoat); 
    
        QE = HC*Aside*(tSC - engine.Twall); 
    end

    % Heat flux out from the convergent part of the nozzle
    Qconv = convergentHeatFlux(At, Ac, engine.thetac, coeff, engine); 
    
    % Temperature of the cooling liquid out from the cooling circuit [K] 
    TF = (engine.fuChamber*QE + engine.fuConv*Qconv)/(mDotF*engine.CPf) + engine.TSf ;        % FUEL
     % + engine.fuConv*Qconv (fattore davanti Qconv)
    TO = (engine.oxChamber*QE + engine.oxConv*Qconv)/(mDotO*engine.CPo) + engine.TSo  ;         % OXIDIZER
    % + engine.oxConv*Qconv (tutto questo termine al numeratore)

    %%% Cost function evaluation
    if TF < engine.TboilF && TO < engine.TboilO
        costVal =  (335 - Ivac/engine.g0)/335;
    else
        costVal = (TF + TO); 
    end

    if Ivac>353*engine.g0
        costVal = abs(costVal)*Ivac/engine.g0;
    end

    
end