function [xc0,xc1,xc2,xc3] = spl100(funf)
global n nn
global deta
global kmat
global alfasp betasp gamasp
global m1 m2
global umat1 lmat1
funspline(1:nn,1:nn-1)=funf;
sm1=(3/deta)*kmat*funspline(2:nn-1,1:nn-1);
sm1(1,1:nn-1)=sm1(1,1:nn-1)+(m1/(2*deta))*funspline(1,1:nn-1);
sm1(n,1:nn-1)=sm1(n,1:nn-1)+(m2/(2*deta))*funspline(nn,1:nn-1);
%funspline_x=zeros(nn,1);
funspline_x(2:nn-1,1:nn-1)=umat1\(lmat1\sm1);
funspline_x(1,1:nn-1)=(1/alfasp)*((betasp/deta)*(funspline(2,1:nn-1)-funspline(1,1:nn-1))... % U_x,0
    +(gamasp/(2*deta)*(funspline(3,1:nn-1)-funspline(1,1:nn-1)))...
    -betasp*funspline_x(2,1:nn-1)-gamasp*funspline_x(3,1:nn-1));
funspline_x(nn,1:nn-1)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1,1:nn-1)-funspline(nn,1:nn-1))... % U_x,N
    -(gamasp/(2*deta))*(funspline(nn-2,1:nn-1)-funspline(nn,1:nn-1))...
    +betasp*funspline_x(nn-1,1:nn-1)+gamasp*funspline_x(nn-2,1:nn-1));

 % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
xc0(1:nn,1:nn-1)=funspline(1:nn,1:nn-1);
xc1(1:nn,1:nn-1)=funspline_x(1:nn,1:nn-1);
xc2(1:nn-1,1:nn-1)=3*(funspline(2:nn,1:nn-1)-funspline(1:nn-1,1:nn-1))/deta^2 ...
            -(2*funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/deta;
xc2(nn,1:nn-1)=0;        
xc3(1:nn-1,1:nn-1)=2*(funspline(1:nn-1,1:nn-1)-funspline(2:nn,1:nn-1))/deta^3 ...
            + (funspline_x(1:nn-1,1:nn-1)+funspline_x(2:nn,1:nn-1))/deta^2;
xc3(nn,1:nn-1)=0;

end

