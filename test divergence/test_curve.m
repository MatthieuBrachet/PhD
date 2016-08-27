function [uu] = test_curve(x,y,z)
    radius=sqrt(x.^2+y.^2+z.^2);
    lambda=atan2(-z,x);
    teta=asin(y./radius);
    
    uu=80*cos(teta).^20;
end

