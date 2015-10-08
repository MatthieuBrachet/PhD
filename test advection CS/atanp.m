function [z] = atanp(y,x)
    z=atan2(y,x);
    if (z<=0) & (z>=-pi)
        z=z+2*pi;
    else
        z=z;
    end
end

