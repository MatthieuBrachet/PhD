function [ up ] = u_test3(tetap)
global u0
xe=.3;
tetab=-pi/6;
tetae=pi/2;
x=xe.*(tetap-tetab).*(tetae-tetab).^(-1);
bx=funb(x);
bxe=funb(xe-x);
up=u0.*bx.*bxe.*exp(4./xe);
end

