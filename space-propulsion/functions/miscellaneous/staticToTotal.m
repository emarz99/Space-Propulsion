function [pT,tT,rhoT] = staticToTotal(gam, Mach,pS, tS,rhoS)
        
    NUM = (1 + 0.5.*(gam-1).*Mach^2); 

    pT = pS.*NUM.^(gam./(gam-1)); 

    tT = tS.*NUM; 

    rhoT = rhoS.*NUM.^(1/(gam-1)); 
end