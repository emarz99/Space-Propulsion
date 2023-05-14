function [mDotOx, mDotFu] = allInj(injectors)

mDotOx = 0 ;
mDotFu = 0 ;

for i=1:lenght(injectors)
    type = injectors(1) ;
    deltaDiameter = injectors(2) ; 
    mDot = singleInj(type, deltaDiameter) ;
    if mDot(1) == 1
        mDotOx = mDotOx + mDot ;
    elseif mDot(1) == 2
        mDotFu = mDotFu + mDot ;
    end
end

end