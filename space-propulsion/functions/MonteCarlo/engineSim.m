function data_engine = engineSim(engine)
% SIMULATE THE PERFORMANCES OF A COMPETE ENGINE, TAKES AS INPUT THE OUTPUT OF designEngine 
%

    %% recall data
    fun = @(x) combustionAnalysis(x, engine);
    
    % Area of a injector [m2] ,  CD of the injector [-], number of injectors [-]
    AinjF = engine.AinjF; CDfu = engine.CDfu; nInjO = engine.nInjO; % Fuel
    AinjO = engine.AinjO; CDox = engine.CDox; nInjF = engine.nInjF; % Oxidizer
    
    rhoF = engine.rhoF; % Fuel density [kg/m3]
    pLossInj = engine.pLossInj;           % Percentage of pressure lost on injection
    
    Ac = engine.Ac;                       % Combustion chamber area [m2]
    At = engine.At;                       % Throat area [m2]
    Ae = At*engine.eps;                       % Exit area [m2]
    Lc = engine.Lc;                       % Chamber length [m]

    mP = engine.mP;                       % total propellant mass [kg]
    

    if length(engine.input.oxidizer.name) == 2
            
        Rgas = 8.3144621; % gas constant [J/kg mol]

        rhoO2 =  (engine.pc*engine.pLossInj)*(0.018)/(engine.input.oxidizer.T(1)*Rgas); 
        rhoH2O = (engine.pc*(1+engine.pLossInj))*(0.032)/(engine.input.oxidizer.T(1)*Rgas); 

        rhoO = (1/3)*rhoO2 + (2/3)*rhoH2O; 

    else
        rhoO = engine.rhoO; 
    end
    %% compute performances 

        % compute the Temperature at the injection plate of the fuel and the
    % combustion chamber pressure
    x0 = engine.pc;
    pc = fmincon(fun, x0, -1, 0); 


    deltaP_inj = pLossInj*pc;                           % pressure drop throug injectino [Pa]

    mDotF = nInjF*AinjF * CDfu*sqrt(2*rhoF*deltaP_inj); % mass flow rate of fuel [kg/s]
    mDotO = nInjO*AinjO * CDox*sqrt(2*rhoO*deltaP_inj); % mass flow rate of oxidizer [kg/s]
    
    mDotP = mDotO + mDotF;                              % propellant mass flow rate [kg/s]
    OF = mDotO/mDotF;                                   % OF ratio
    tb = mP/mDotP;                                      % Burning time

    %%% CEA input
    engine.input.OF = OF; 
    engine.input.pc = pc*(1e-5);
    engine.input.supar = engine.eps; 
    engine.eps = Ae/At;
    coeff = createDataCEA(engine.input, engine.ceaPath); 

    % Recall coefficients
    p = coeff(1, :);    
    Temp = coeff(2, :); 
    rho = coeff(3, :);   
    gam = coeff(4, :); 
    vSon = coeff(5, :);
    Mach = coeff(6, :); 
    
    % Compute nozzle parameters
    [lambda, lNoz, thetaI, thetaE] = nozzleDesign(At, engine); 
    
    Thrust = lambda*mDotP*Mach(end)*vSon(end) + p(end)*Ae;      % Engine thrust [N]

    Isp = Thrust/(mDotP*engine.g0);                             % Engine specific impulse at vac [s]
    
    deltaV = Isp*engine.g0*log(1 + mP/engine.ms);                   % Possible deltaV
    
    [delta_pTankOx, delta_pTankFu] = TotalPressureDrop(true, true, mDotO, mDotF, engine); 
    pTankFu = delta_pTankFu + pc;
    pTankOx = delta_pTankOx + pc;
    
    % Heat transfer through convergent
    [Qconv,qConv, qConvSpec] = convergentHeatFlux(At, Ac, engine.thetac, coeff, engine);
    
    % Compute flow velocity in the chamber
    f = @(x) (mDotP/(Ac*rho(1)))*(1 + 0.5*(gam(1)-1)*(x/vSon(1))^2)^(1/(gam(1)-1)) - x; 
    vC = fzero(f, 0);
    machC = vC/vSon(1); 

    % obtain static temperature in the chamber
    [~, tSC, ~] = totalToStatic(gam(1), machC, p(1), Temp(1), rho(1));

    Rc = sqrt(Ac/pi);       % Chamber radius [m]
    Aside = 2*Rc*pi*Lc;     % Chamber side area [m2]
    
    hC = hConvection(2*Rc, vC, coeff(:, 1), true); 

    % Heat flux out from the combustion chamber [W]
    %QE = hC*Aside*(tSC(1) - engine.Twall); 
    
    if isnan(engine.kCoat)
        QE = hC*Aside*(tSC(1) - engine.Twall); 
    else
        xIs = ((tSC - engine.Twall)/(tSC - engine.Tcoat) - 1) * engine.kCoat/hC;
    
        HC = 1/(1/hC + xIs/engine.kCoat); 
    
        QE = HC*Aside*(tSC - engine.Twall); 
    end  
    
    qSpec = QE/Aside; 

    qConv = [QE, qConv]; 
    qConvSpec = [qSpec, qConvSpec];
    
    
    % Temperature of the cooling liquid out from the cooling circuit [K] 
    TF = (engine.fuChamber*QE + engine.fuConv*Qconv)/(mDotF*engine.CPf) + engine.TSf;  % FUEL
     
    TO = (engine.oxChamber*QE + engine.oxConv*Qconv )/(mDotO*engine.CPo) + engine.TSo; % OXIDIZER


    % Compute the injection velocity [m/s]
    vOinj = CDox * sqrt(2*pc*pLossInj/rhoO);    % Oxidizer
    vFinj = CDfu * sqrt(2*pc*pLossInj/rhoF);    % Fuel

    %% Save data
    data_engine.Thrust = Thrust; 
    data_engine.Isp = Isp; 
    data_engine.deltaV = deltaV; 
    data_engine.tb = tb; 

    data_engine.mDotP = mDotP; 
    data_engine.mDotO = mDotO; 
    data_engine.mDotF = mDotF; 
    
    data_engine.mP = mP;
    data_engine.mF = engine.mF; 
    data_engine.mO = engine.mO; 
    
    data_engine.OF = OF; 
    data_engine.pc = pc; 
    
    data_engine.coeff = coeff; 
    
    data_engine.qConv = qConv; 
    data_engine.qConvSpec = qConvSpec; 
    
    data_engine.TF = TF;
    data_engine.TO = TO; 

    data_engine.Tflame = Temp(1); 
    data_engine.nozzle.thetaI = thetaI; 
    data_engine.nozzle.thetaE = thetaE; 
    data_engine.nozzle.Ln = lNoz; 
    data_engine.nozzle.lambda = lambda; 
    data_engine.pTankFu = pTankFu; 
    data_engine.pTankOx = pTankOx;
    data_engine.vOinj = vOinj; 
    data_engine.vFinj = vFinj;

end
