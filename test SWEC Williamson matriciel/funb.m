function [ bx ] = funb( x )
for i=1:size(x,1)
    for j=1:size(x,2)
        if 0<x(i,j)
            bx(i,j)=exp(-1./x(i,j));
        else
            bx(i,j)=0;
        end
    end
end
end

