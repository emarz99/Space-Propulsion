function [qTot,q,qSpec] = convergentHeatFlux(At, Ac, thetaC, coeff,engine)
% Compute the heat flux out from the convergent part of the nozzle
%
% qTot: total heat flux [W]
% q: heat flux at each station [W]
% qSpec: specific heat flux at each station [W/m2]


    SlatFun = @(r1, r2, k) pi.*(r1+r2).*sqrt(k.^2 + (r1 - r2).^2 ); 
    getCoeff = @(i) coeff(i, [3:2+engine.nSub, 2]);


    A = [Ac engine.input.subar.*At At];  

    r = sqrt(A./pi);
    
    h = [0, (r(1)-r(2:end))./tan(thetaC)]; 

    Slat= SlatFun(r(1:end-1), r(2:end), diff(h)); 
    
    Mach = getCoeff(6);  
    Temp = getCoeff(2);
    vSon = getCoeff(5);

    vC = vSon.*Mach;   

    
    hC = hConvection(2*r(1:end-1), vC, coeff(:, [3:2+engine.nSub, 2]), false); 

    q = hC.*Slat.*(Temp - engine.Twall); 
    qSpec = q./Slat; 

    qTot = sum(q); 


end