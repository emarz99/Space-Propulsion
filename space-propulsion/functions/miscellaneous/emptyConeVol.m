function [volume] = emptyConeVol(radiusMax, radiusMin, lenght, thickness)

% truncated cone with both bases opened

extVolume = volConeCut(lenght, radiusMax, radiusMin) ;
intVolume = volConeCut(lenght, RadiusMax-thickness, radiusMin-thickness) ; 

volume = extVolume - intVolume  ;

end