function [velocity] = radius2velocity(massFlowRate, Rho, Radius)

Area = pi*Radius^2 ; 

velocity = massFlowRate/(Rho*Area); 

end