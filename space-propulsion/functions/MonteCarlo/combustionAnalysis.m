function delta = combustionAnalysis(x, engine)
    
% find the correct ombination of Tinj Fuel and combustion chamber pressure
% with the given constraints

%  x(1) -> fuel temperature at injection
%  x(2) -> combustion pressure

    deltaP_inj = engine.pLossInj*x;  % Injection pressure loss
    engine.eps = engine.Ae/engine.At;   % Expansion ratio
    
    if length(engine.input.oxidizer.name) == 2
            
        Rgas = 8.3144621; % gas constant [J/kg mol]

        rhoO2 =  (engine.pc*engine.pLossInj)*(0.018)/(engine.input.oxidizer.T(1)*Rgas); 
        rhoH2O = (engine.pc*(1+engine.pLossInj))*(0.032)/(engine.input.oxidizer.T(1)*Rgas); 

        rhoO = (1/3)*rhoO2 + (2/3)*rhoH2O; 

    else
        rhoO = engine.rhoO; 
    end


    % Mass flow rate [kg/s]
    mDotF = engine.nInjF*engine.AinjF * engine.CDfu*sqrt(2*engine.rhoF*deltaP_inj); % Fuel
    mDotO = engine.nInjO*engine.AinjO * engine.CDox*sqrt(2*rhoO*deltaP_inj); % Oxidizer
    mDotP = mDotO + mDotF;                                                          % total propellant
    
    %%% CEA input
    engine.input.pc = x*(1e-5);
    engine.input.supar = engine.eps; 
    engine.input.OF = mDotO/mDotF;

    %%% Recall coefficients
    coeff = createDataCEA(engine.input, engine.ceaPath); 
    p = coeff(1, :);    
    Temp = coeff(2, :); 
    rho = coeff(3, :);   
    gam = coeff(4, :); 
    vSon = coeff(5, :);
    cStar = coeff(9, :); 

    

    Ac = engine.Ac;                     % chamber area [m2]
    At = engine.At;                     % throat  area [m2]
    Rc = sqrt(Ac/pi);                   % chamber radius [m]
    Lc = engine.Lc;                     % chamber length [m]
    

    Aside = 2*Rc*pi*Lc; 
    
    % compute flow velocity in the chamber to satisfy the continuity
    % equation
    f = @(x) (mDotP/(Ac*rho(1)))*(1 + 0.5*(gam(1)-1)*(x/vSon(1))^2)^(1/(gam(1)-1)) - x; 
    
    vC = fzero(f, 0);                           % velocity inside the chamber [m/s]
    machC = vC/vSon(1);                         % mach inside the chamber     

    % obtain static temperature in the chamber
    [~, tSC, ~] = totalToStatic(gam(1), machC, p(1), Temp(1), rho(1));
    
    hC = hConvection(2*Rc, vC, coeff(:, 1), true);    % compute convection heat flux coefficient [W/(m2 K)]
    
    % Heat flux out from the combustion chamber [W]
    QE = hC*Aside*(tSC(1) - engine.Twall); 
    
    % Heat flux out from the convergent part of the nozzle
    Qconv = convergentHeatFlux(At, Ac, engine.thetac, coeff, engine); 
    
    %%% Chamber pressure
    pc = mDotP/At * cStar(end); 

    
%    delta(1) = abs((QE+Qconv)/(mDotF*engine.CPf) + engine.TS - x(1))/x(1);
    delta= abs((pc - x)/x); 

end