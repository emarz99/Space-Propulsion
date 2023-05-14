function [AreaPipe] = AreaP(MassFlowRate, FlowVelocity, Rho)

AreaPipe = MassFlowRate./(FlowVelocity.*Rho) ;

end