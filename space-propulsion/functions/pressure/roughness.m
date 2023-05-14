function epsilon = roughness(varargin)
%ROUGHNESS  Compute the relative roughness coefficient of a pipe from
%   values of the friction factor and Reynolds number for different
%   operating points.
%-----------------------------------------------------------
%   Syntax:
%      epsilon = roughness(Re,f)
%      epsilon = roughness(Hin,Hout,Q,nu,L,D,g)
%-----------------------------------------------------------
%           Re : Reynolds number
%            f : Darcy-Weisbach friction factor
%          Hin : Piezometric head at pipe inlet [m]
%         Hout : Piezometric head at pipe outlet [m]
%            Q : Flow rate [m^3/s]
%           nu : Kinematic viscosity [m^2/s]
%            L : Pipe length [m] - scalar value
%            D : Pipe diameter [m] - scalar value
%            g : Gravity acceleration [m/s^2] - scalar value
%      epsilon : Relative roughness coefficient
%-----------------------------------------------------------
%   Example:
%      Re = [47525, 74725, 99490, 123013];
%      f = [0.022786, 0.021086, 0.020241, 0.019698]; 
%      epsilon = roughness(Re,f)
%-----------------------------------------------------------
%   Author:
%      Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx
%-----------------------------------------------------------
%   Cite as:
%      Santos-Ruiz, Ildeberto. (2018, November 9).
%      Friction and Roughness. Zenodo.
%      http://doi.org/10.5281/zenodo.1481562
%-----------------------------------------------------------

switch nargin
    case 2
        % roughness(Re,f)
        Re = varargin{1};
        f = varargin{2};
    case 7
        % roughness(Hin,Hout,Q,nu,L,D,g)
        Hin = varargin{1};
        Hout = varargin{2};
        Q = varargin{3};
        nu = varargin{4};
        L = varargin{5};
        D = varargin{6};
        g = varargin{7};
        Re = 4*Q./(pi*nu*D);
        f = g*pi^2*D^5*(Hin-Hout)./(8*L*Q.^2);
    otherwise
        error('Only 2 or 7 input arguments are accepted!')
        return
end

colebrook = @(f,Re,epsilon) 1./sqrt(f)+2*log10(epsilon/3.7+...
    2.51./(Re.*sqrt(f)));
objfunc = @(epsilon) colebrook(f(:),Re(:),epsilon);
options = optimoptions('lsqnonlin','Display','none');
epsilon = lsqnonlin(objfunc,eps,0,1,options);

end