function [Re] = Reynolds(Rho, FlowVelocity, Diameter, Mu)

Re = (Rho.*FlowVelocity.*Diameter)./Mu ;

end