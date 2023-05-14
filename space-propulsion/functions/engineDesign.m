
function engineFinal = engineDesign(xOpt, engine)    
% USES THE OUTPUT OF THE OPTIMIZER TO DEFINE THE COMPLETE ENGINE AND SAVES IN ENGINEFINALE    
%

    %% Recall data from optimization
    engine.OF = xOpt(1);

    Tobj = xOpt(2); 
    
    %% CEA run

    % create CEA input
    engine.input.OF = xOpt(1); 
    engine.input.subar = linspace(engine.CR, 1, engine.nSub+2);
    engine.input.supar = linspace(1, engine.eps, engine.nSup+1); 

    engine.input.subar = engine.input.subar(2:end-1); 
    engine.input.supar = engine.input.supar(2:end); 
    
    % recall output
    coeff = createDataCEA(engine.input, engine.ceaPath); 
    p = coeff(1, :);    
    Temp = coeff(2, :); 
    rho = coeff(3, :);   
    gam = coeff(4, :); 
    vSon = coeff(5, :);
    Mach = coeff(6, :); 
    cStar = coeff(9, :); 
    cT = coeff(10, :);

    

    %% Computation -------------------------------------------------------
    
    %%% Same as engineOpt

    % function to compute the throat area to obtain the given thrust
    AtFun = @(x) nozzleDesign(x, engine)*((engine.pc*x)/cStar(end))*Mach(end)*vSon(end) + p(end)*x*engine.eps - Tobj; 
    At0 = Tobj/(cT(end)*engine.pc);    % initial guesss throat area [m^2]
       
    [At, res] = fzero(AtFun, At0); 

    mDotP = (engine.pc*At)/cStar(end);          % propellant mass flow rate [kg/s]
    mDotF = mDotP/(engine.OF + 1);              % fuel mass flow rate [kg/s]
    mDotO = mDotF*engine.OF;                    % oxidizer mass flow rate [kg/s]

    Isp = Tobj/(mDotP*engine.g0);                           % specific impulse [m/s]

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
    %QE = hC*Aside*(tSC(1) - engine.Twall); 
    
    if isnan(engine.kCoat)
        QE = hC*Aside*(tSC(1) - engine.Twall); 
    else
        xIs = ((tSC - engine.Twall)/(tSC - engine.Tcoat) - 1) * engine.kCoat/hC;
    
        HC = 1/(1/hC + xIs/engine.kCoat); 
    
        QE = HC*Aside*(tSC - engine.Twall); 
    end

    % Heat flux out from the convergent part of the nozzle
    [Qconv , qVecConv, qSpecConv] = convergentHeatFlux(At, Ac, engine.thetac, coeff, engine); 
    
    % Temperature of the cooling liquid out from the cooling circuit [K] 
    TF = (engine.fuChamber*QE + engine.fuConv*Qconv)/(mDotF*engine.CPf) + engine.TSf;  % FUEL
     
    TO = (engine.oxChamber*QE + engine.oxConv*Qconv )/(mDotO*engine.CPo) + engine.TSo; % OXIDIZER
    
    %% Injection plate
    
    if length(engine.input.oxidizer.name) == 2
            
        Rgas = 8.3144621; % gas constant [J/kg mol]

        rhoO2 =  (engine.pc*engine.pLossInj)*(0.018)/(engine.input.oxidizer.T(1)*Rgas); 
        rhoH2O = (engine.pc*(1+engine.pLossInj))*(0.032)/(engine.input.oxidizer.T(1)*Rgas); 

        rhoO = (1/3)*rhoO2 + (2/3)*rhoH2O; 

    else
        rhoO = engine.rhoO; 
    end


    AinjF_tot = mDotF/(engine.CDfu*sqrt(2*engine.rhoF*engine.pLossInj*engine.pc)); 
    AinjO_tot = mDotO/(engine.CDox*sqrt(2*rhoO*engine.pLossInj*engine.pc));

    [Ainj, nInj] = injectorsDesign(AinjF_tot, AinjO_tot, engine); 
    
    vOinj = engine.CDox * sqrt(2*engine.pc*engine.pLossInj/engine.rhoO);    % Oxidizer injection velocity
    vFinj = engine.CDfu * sqrt(2*engine.pc*engine.pLossInj/engine.rhoF);    % Fuel injection velocity
    %% Tanks

    [deltaP_ox, deltaP_fu, lMaxOx, lMaxFu] = TotalPressureDrop(true, true, mDotO, mDotF, engine); 
    pTankO = engine.pc + deltaP_ox;  % [Pa]
    pTankF = engine.pc + deltaP_fu;  % [Pa]
    
    %% Masses

    mP = engine.ms*( exp(engine.deltaV/(Isp*engine.g0)) - 1 );  % Total propellant mass [kg]
    mF = mP/(1 + engine.OF);                                    % Fuel mass [kg]
    mO = engine.OF*mF;                                          % Oxidizer mass [kg]

    tb = mP/mDotP;                                              % Burning time [s]
    
    volF = mF/engine.rhoF; 
    volO = mO/engine.rhoO; 

    [vH, mH] = heliumTankDesign(pTankF, pTankO, volF, volO, engine);
    %% Save data
    
    engineFinal = engine; 
    
    engineFinal.Thrust = Tobj; 
    engineFinal.Isp = Isp; 
    engineFinal.At = At; 
    engineFinal.tb = tb; 
    engineFinal.Ae = At*engine.eps; 
    engineFinal.Tflame = Temp(1); 
    engineFinal.pt = p(2); 
    engineFinal.pe = p(end); 

    engineFinal.Ac = Ac; 
    engineFinal.Lc = Lc; 
    
    engineFinal.mP = mP; 
    engineFinal.mF = mF;
    engineFinal.mO = mO;

    engineFinal.mDotP = mDotP; 
    engineFinal.mDotF = mDotF; 
    engineFinal.mDotO = mDotO; 
    
    engineFinal.TF = TF; 
    engineFinal.TO = TO;
    if not(isnan(engine.kCoat))
        engineFinal.xIS = xIs; 
    end

    engineFinal.convectHeatFlux = [QE, qVecConv]; 
    engineFinal.convectSpecHeatFlux = [QE/Aside, qSpecConv]; 

    engineFinal.AinjF = Ainj(1); 
    engineFinal.AinjO = Ainj(2); 
    
    engineFinal.nInjF = nInj(1); 
    engineFinal.nInjO = nInj(2); 

    engineFinal.vInjF = vFinj; 
    engineFinal.vInjO = vOinj; 

    engineFinal.tankPressureF = pTankF; 
    engineFinal.tankPressureO = pTankO; 
    engineFinal.tankVolumeF = volF; 
    engineFinal.tankVolumeO = volO; 
    engineFinal.tankVolumeHe = vH; 
    engineFinal.mH = mH; 

    engineFinal.maxFuelPipeLength = lMaxFu; 
    engineFinal.maxOxPipeLength = lMaxOx; 
    engineFinal.flowVelocityOx = mDotO/(engine.rhoO*engine.pipeDiamO^2*pi/4); 
    engineFinal.flowVelocityFu = mDotF/(engine.rhoF*engine.pipeDiamF^2*pi/4); 
end

