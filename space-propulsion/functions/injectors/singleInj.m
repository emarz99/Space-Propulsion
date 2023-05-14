function [mDot, type] = singleInj(type, deltaDiameter)

% type account for fu or ox
% type == 1 ox type == 2 fu
cdO = 1 ;
cdF = 1 ;

if type == 1
    rho = engine.rhoO ;
    diameter = diameterO ;
elseif type == 2
    rho = engine.rhoF ;
    diameter = diameterF ;
end

deltaP = engine.pc*0.3 ;
DeltaArea = ((diameter + deltaDiameter)^2) ;

mDot = rho*Cd*pi*DeltaArea*((2*deltaP/rho)^0.5) ;

end