% SHOW THE RATE OF CONVERGENCE
% OF THE INTEGRATION FORMULA
% FOR A CHOSEN (nhs,mhs) SPH

clear all;
global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global dxi deta dga;

%%%%%%%%%%%%%%%%%%%

% Number of experiments
k=5;

nhs=100;
mhs=4;
err_i=zeros(k,1);
ii_a=zeros(k,1);
kN=1;
for N=2.^(2:k+1)
    make_cs_grid(N); % APPEL DU MODULE 74
    %eps_w=zeros(nn,nn);
    [eps_w,opt_val]=eps_weights(20);
    % FUNCTION TO INTERPOLATE - WILL BE GIVEN IN THE ARGUMENTS.
    funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
    funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
    funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
    funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
    funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
    funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);
    % APPROXIMATE QUADRATURE ON THE SPHERE
    [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
        int_weights(dxi*deta*(dga+eps_w),funfI,funfII,funfIII,funfIV,funfV,funfVI);
    ii_a(kN)=nrmg;
    ii_ex=0; % EXACT INTEGRAL OF SPHERICAL HARMONICS FOR N.=1, M>=1;
    err_i(kN)=abs(ii_a(kN)-ii_ex);
    kN=kN+1;
end
err_i
%plot(log(2.^(2:(k+1))),-log(err_i));