function [Ymn] = spharm(L,M,X,Y,Z)
% Build the spherical harmonic on the point of (X,Y,Z) for 
% L         - Spherical harmonic degree, [1x1]
% M         - Spherical harmonic order,  [1x1]
% Define constants (REQUIRED THAT L(DEGREE)>=M(ORDER)).

% Author : Matthieu Brachet (IECL)
% DBE 2016/10/19

[LAMBDA,TETA,R]=cart2sph(X,Y,Z);
Lmn=legendre(L,cos(TETA-pi/2));

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Lmn=squeeze(Lmn(M+1,:,:));
Ymn=C*Lmn.*exp(1i*M*LAMBDA);
end