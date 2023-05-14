function [lambda, L, thetaI, thetaE] = nozzleDesign(At, engine)
    
    getTheta = @(data) interp1(data(:, 1), data(:, 3), engine.eps)*pi/180;


    % consider max 0.8L (wrt nozzle at 15deg)
    Rt = sqrt(At/pi);
    L = 0.8*(sqrt(engine.eps) - 1)*Rt/tand(15); 
    
    
    alpha = atan(1/0.8 * tan(15*pi/180));


    thetaI = getTheta(engine.RAOinitial);
    thetaE = getTheta(engine.RAOfinal);
    
    lambda = 0.5*(1 + cos( (alpha + thetaE)/2 ) );
    

end