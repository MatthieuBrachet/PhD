function [z] = solver_impli75(bb,zeta,ddt,tau,x,y,z)
%% *** RSS SOLVER ****************************************************
% solve :
% z + ddt * (2*omega*cos(teta) + zeta)*k \wedge z = bb
%on each point x,y,z.
global radius omega
[n1,n2,~]=size(bb);
[~, teta,~]=cart2sph(x,y,z);
for i=1:n1
    for j=1:n2
        b(1:3)=bb(i,j,1:3);
        kx=x(i,j)./radius;
        ky=y(i,j)./radius;
        kz=z(i,j)./radius;
        PDT=[0 -kz ky;kz 0 -kx; -ky kx 0 ];
        A=speye(3,3)+tau*ddt*(2*omega*sin(teta(i,j))+zeta(i,j))*PDT;
        z(i,j,1:3)=A\b';
    end
end

end
