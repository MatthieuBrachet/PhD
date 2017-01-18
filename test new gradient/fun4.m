function [ fun ] = fun4( x,y,z )
global radius
% fun=x.^3+y.^2+z;
% r=radius^3+radius^2+radius;
% fun=fun./r;
nhs=100; mhs=2;
fun=sph( nhs,mhs,x,y,z );
fun=imag(fun);
end