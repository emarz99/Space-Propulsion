function [volume] = volConeCut(lenght, RadiusMax, RadiusMin)

volume = 1/3*pi*lenght*(RadiusMax^2 + RadiusMin^2 + RadiusMax*RadiusMin )  ;

end