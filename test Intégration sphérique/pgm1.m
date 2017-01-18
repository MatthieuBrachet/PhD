%%%%% SHOW THE INTEGRATION ERRORS  
%%%%% FOR SEVERAL n (LINES) AND m (COLUMNS)

clear all;
global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

%%%%%%%%%%%%%%%%%%%

% TEST FOR QUADRATURE FORMULA IN NRM74: EMPIRICAL QUADRATURE
% WITH COEFFICIENTS 1/2 (EDGE) AND 1/3 (CORNER).
mod74; % APPEL MODULE "PROBLEME"
%
% CASE 1; QUADRATURE OF SPHERICAL HARMONICS
% --------------------------------------------
% m1 and m2 multiples of 4
m1=4;
m2=12;
% n1 and n2 multiples of 2
n1=12;
n2=200;
err_i=zeros((n2-n1)/2+1,(m2-m1)/4+1);
ii_a=zeros((n2-n1)/2+1,(m2-m1)/4+1);
kn=1;
for nhs=n1:2:n2
    km=1;
    for mhs=m1:4:m2
        % FUNCTION TO INTERPOLATE - WILL BE GIVEN IN THE ARGUMENTS.
        funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
        funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
        funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
        funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
        funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
        funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);
        % APPROXIMATE QUADRATURE ON THE SPHERE
        str='int';
        [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
            nrm74(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);
        ii_a(kn,km)=nrmg;
        ii_ex=0; % EXACT INTEGRAL OF SPHERICAL HARMONICS FOR N.=1, M>=1;
        err_i(kn,km)=abs(ii_a(kn,km)-ii_ex);
        km=km+1;
    end
    kn=kn+1;
end
err_i