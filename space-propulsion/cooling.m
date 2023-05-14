% Coolant-side and Gas-side Properties

% Coolant-side and Gas-side Properties 1 - Propane (Liquid) - Non Toxic
% These data need to be double-checked
rho_f =  ;  % density of fuell at 15°C and 1 atm [kg/m^3]
mu_f =  ;      % viscosity of fuell
k_f =  ;  % heat conductivity of fuell [W/m K]
c_f =   ;  % specific heat of fuell 

% Coolant-side and Gas-side Properties 2 - MMH - Toxic

% rho_f = ;  % density of fuell 
% mu_f =  ;  % viscosity of fuell
% k_f =   ;  % heat conductivity of fuell
% c_f =   ;  % specific heat of fuell 

% Material for Combustion Chamber Properties  

% Material Properties 1 - Niobium Alloy C-103 - Non Toxic 

t_max = 2350 ;       % melting point [°C]
rho = 8.19*(10^3) ;  % density [kg/m^3]
k = 37.4 ;           % thermal conductivity at 1600°C [W/m K]
% k = 42.4 ;         % thermal conductivity at 800°C [W/m K]
c = 340 ;            % specific heat [J/kg k]

% Material Properties 2 - Inconel Alloy 718 - Toxic
% minimum hole diameter 1mm
% t_max = 2300 ;       % melting point [°C]
% rho = 8.19*(10^3) ;  % density [kg/m^3]
% k = 11.4 ;           % thermal conductivity [W/m K]
% c = 435 ;            % specific heat [J/kg k]

% Combustion Chamber Specifications

t = 700*10^(-6) ;         % thickness [m] (first guess)
% r = ;                     % radius of cooling channels [m] (first guess)
d = 2*r  ;                % diameter of cooling channels [m] (first guess)
L = (4*pi*r^2)/2*pi*r  ;  % hydraulic diameter of cooling channels [m] (first guess)

% Computing Important Parameters

pr = (cp_f*mu_f)/k_f ;              % Prantl Number
% re = (rho_f* *L)/mu_f ;             % Reynolds Number
n_h = 0.4 ;                         % Nusselt Pr Coefficient for hot side
n_c = 0.3 ;                         % Nusselt Pr Coefficient for cool side
nu_h = 0.023*(re^(4/5))*(pr^n_h) ;  % Nusselt Number for heated side
nu_c = 0.023*(re^(4/5))*(pr^n_c) ;  % Nusselt Number for cooled side

h_h = (nu_h*k_f)/L ; % Hot-Side convective heat transfer
h_c = (nu_c*k_f)/L ; % Coolant-Side convective heat transfer

h = ((1/h_h)+(t/k)+(1/h_c))^(-1) ;


% check if Dittus-Boelter Equation is applicable
if pr > 0.6 || pr<160 || re>10000 
    disp('Dittus-Boelter Equation is applicable');
else
    disp('Dittus-Boelter Equation is not applicable')
end


