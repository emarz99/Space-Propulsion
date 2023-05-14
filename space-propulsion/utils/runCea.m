

ceaPath = '../CEA/'; 
addpath(ceaPath);

x2=CEA('problem','rocket','frozen','nfz',2,'o/f',2.5,'p(bar)',10,'subsonic(ae/at)',3,2,'supsonic(ae/at)',[1,40],'reactants','fuel','N2H4','t(k)',300.0,'oxid','N2O4','t(k)',300,'output','transport','end','screen');
