function [volume] = emptyCilinderVol(shape, radius, lenght, thickness)
% 1 is tank
% tank is a complitely closed cilinder 
% 2 is pipe
% pipe is a cilinder with base areas open
% 3 is combustion chamber 
% combustion chamber is a cilinder with one base area and one base area
% closed
% volume of cilindrical tank given ext radius, thickness, lenght

extVolume = cilinderVol(radius, lenght) ;
intVolume = cilinderVol(radius-thickness, lenght) ; 

if shape == 1 
    baseVolume = 2*cilinderVol(radius-thickness, thickness) ;
elseif shape == 2 
    baseVolume = 0 ;
elseif shape == 3
    baseVolume = cilinderVol(radius-thickness, thickness) ;
end

volume = extVolume - intVolume + (baseVolume) ;

end