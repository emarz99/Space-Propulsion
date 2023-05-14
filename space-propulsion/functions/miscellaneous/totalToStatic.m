function [pS, tS, rhoS] = totalToStatic(gam, Mach,pT, tT, rhoT)

    DEN = (1 + 0.5.*(gam-1).*Mach.^2); 

    pS = pT./DEN.^(gam./(gam-1)); 

    tS = tT./DEN; 
    
    rhoS = rhoT./(DEN.^(1/gam-1));
end

