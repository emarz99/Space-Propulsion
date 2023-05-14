function [vHe, mHe] = heliumTankDesign(pTankFu, pTankOx, vTankFu, vTankOx, engine)
    gam = 5/3;     % monoatomic gas
    R = 8.3144621; % gas constant [J/kg mol]
    
    pDim = max(pTankFu, pTankOx);
    rhoH = (engine.pHelium*engine.MMhelium)/(R*engine.TShelum); 

    NUM = (vTankOx+ vTankFu); 
    DEN = (engine.pHelium/pDim)^(1/gam) - 1; 

    vHe = NUM/DEN; 

    mHe = vHe*rhoH; 


end

