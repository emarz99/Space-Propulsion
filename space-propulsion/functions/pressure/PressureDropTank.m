function [DeltaPressure] = PressureDropTank(Rho, FlowVelocity)

DeltaPressure = 0.5.*Rho.*(FlowVelocity).^2 ; % Pressure Drop for Tank Exit

end