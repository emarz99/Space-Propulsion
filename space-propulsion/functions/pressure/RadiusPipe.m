function [Radius] = RadiusPipe(MassFlowRate, FlowVelocity, Rho)

AreaPipe = MassFlowRate./(FlowVelocity.*Rho) ;

Radius = (AreaPipe/pi).^(0.5) ;

end