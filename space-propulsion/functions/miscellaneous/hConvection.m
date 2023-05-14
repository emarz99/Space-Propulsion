function h = hConvection(L, vel, coeff, isTotal)
% Compute che convection heat transfer coefficient


    p = coeff(1, :);    
    Temp = coeff(2, :); 
    rho = coeff(3, :);   
    gam = coeff(4, :); 
    vSon = coeff(5, :);
    kCond = coeff(7, :);
    Pr = coeff(8, :);
    visc = coeff(13, :); 
    
    % convert total quantities to static ones
    if isTotal
        machC = vel./vSon; 
        [~, ~, rho] = totalToStatic(gam, machC, p, Temp, rho); 
    end
    Re = (rho.*vel.*L)./(visc);
    
% %   ***** Gnielinski ***** 
%     K = 0.002 ; % Absolute Roughness Coefficient for Alluminum [mm]
%     Epsilon = K./(L) ; % Relative Roughness Coefficient for Alluminum Oxydizer [mm]
%     Darcy = friction(Re, Epsilon, 'iterative') ;
%     Nu = ((Darcy./8).*(Re - 1000).*Pr)./(1 + 12.7*(Darcy./8).^(0.5).*(Pr.^(2/3) - 1) ); 
    
%  ***** Dittus Boelter
   Nu = (0.0265.*Re.^0.8.*Pr.^0.3);  
    
    h = (kCond./L).*Nu; 
end

