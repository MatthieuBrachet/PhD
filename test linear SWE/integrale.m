function [ gp ] = integrale( phi )
global g radius OMEGA

syms x
for i=1:size(phi,1)
    for j=1:size(phi,2)
        pp=phi(i,j);
        gp(i,j)=g*h0-double(int(radius*u_test(x)*(2*OMEGA*sin(x)+tan(x)./radius.*u_test(x)),0,pp));
    end
end


end

