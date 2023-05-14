function [DeltaPressure] = PressureDropFeedLine(Darcy, Rho, FlowVelocity, Length, Diameter)

DeltaPressure = 0.5.*Darcy.*Rho.*((FlowVelocity).^2).*Length.*Diameter ; % Pressure Drop for Feed Line

end