function f = friction(Re,epsilon,method)
%FRICTION  Computes the Darcy-Weisbach friction factor for
%   turbulent flow in pipes using different methods.
%----------------------------------------------------------
%   Syntax:
%      f = friction(Re,epsilon,method)
%----------------------------------------------------------
%           Re : Reynolds number.
%      epsilon : Relative roughness coefficient, k/D.
%       method : Algorithm or method to use.
%                Available methods are:
%                'Swamee-Jain', 'Haaland', 'Serghides',
%                'LambertW' and 'iterative' (default).
%            f : Darcy-Weisbach friction factor.
%----------------------------------------------------------
%   Examples:
%      f = friction(1e5,1e-4)
%      f = friction(1e5,1e-4,'Swamee-Jain')
%      [Re,epsilon] = meshgrid(1e4:1e4:1e5,1e-4:1e-4:1e-3);
%      f = friction(Re,epsilon,'Serghides');
%      surf(Re,epsilon,f)
%----------------------------------------------------------
%   Author:
%      Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx
%----------------------------------------------------------
%   Cite as:
%      Santos-Ruiz, Ildeberto. (2018, November 9).
%      Friction and Roughness. Zenodo.
%      http://doi.org/10.5281/zenodo.1481562
%---------------------------------------------------------
if nargin < 3
    method = 'iterative';
end
switch method
    case 'Swamee-Jain'
        f = SwameeJain(Re,epsilon);
    case 'Haaland'
        f = Haaland(Re,epsilon);
    case 'Serghides'
        f = Serghides(Re,epsilon);
    case 'LambertW'
        f = LambertW(Re,epsilon);
    case 'iterative'
        f = zeros(size(Re));
        for i = 1:size(Re,1)
            for j = 1:size(epsilon,2)
                f(i,j) = ColebrookWhite(Re(i,j),epsilon(i,j));
            end
        end
end
    function f = SwameeJain(Re,epsilon)
        f = 0.25./log10(5.74./Re.^0.9+epsilon/3.7).^2;
    end
    function f = Haaland(Re,epsilon)
        f = 1./(-1.8*log10(6.9./Re+(epsilon/3.7).^1.11)).^2;
    end
    function f = Serghides(Re,epsilon)
        A = -2*log10(epsilon/3.7+12./Re);
        B = -2*log10(epsilon/3.7+2.51*A./Re);
        C = -2*log10(epsilon/3.7+2.51*B./Re);
        f = 1./(A-(B-A).^2./(C-2*B+A)).^2;
    end
    function f = LambertW(Re,epsilon)
        a = 2.51./Re;
        b = epsilon/3.7;
        f = 1./(2*lambertw(0,log(10)./(2*a).*10.^(b./(2*a)))/log(10)-b./a).^2;
    end
    function f = ColebrookWhite(Re,epsilon)
        fun = @(f) 1./sqrt(f)+2*log10(epsilon/3.7+2.51./(Re.*sqrt(f)));
        f = fzero(fun,[eps,1]);
    end
end