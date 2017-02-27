clc; clear all; close all;
% fun1 is spherical harmonic arounf (Oz).
global n nn
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=32;
make_cs_grid(N);
axes='oz';

% Computation of weights with basic formula
weights=dxi*deta*dga;

% computation of weights with B. Portenelle's formula
% disp('weights...')
% nhs_max=100;
% k=0;
% [A,err_i]=compute_A_sym(nhs_max);
% eps_w=solve_weights(A,err_i,k,1,1);
% weights=dxi*deta*(dga+eps_w);




nhs1=11; mhs1=5;
fun1fI=conj(sph(nhs1,mhs1,x_fI,y_fI,z_fI));
fun1fII=conj(sph(nhs1,mhs1,x_fII,y_fII,z_fII));
fun1fIII=conj(sph(nhs1,mhs1,x_fIII,y_fIII,z_fIII));
fun1fIV=conj(sph(nhs1,mhs1,x_fIV,y_fIV,z_fIV));
fun1fV=conj(sph(nhs1,mhs1,x_fV,y_fV,z_fV));
fun1fVI=conj(sph(nhs1,mhs1,x_fVI,y_fVI,z_fVI));

%nhs2=nhs1;
nhs_max=25;%2*nhs1;
for i=1:nhs_max
    nhs2=i-1;
    nhs_test=-nhs2:1:nhs2;
    for j=1:length(nhs_test)
        mhs2=nhs_test(j);
        clc; mhs2

        if strcmp(axes,'oz')==1
            fun2fI=sph(nhs2,mhs2,x_fI,y_fI,z_fI);
            fun2fII=sph(nhs2,mhs2,x_fII,y_fII,z_fII);
            fun2fIII=sph(nhs2,mhs2,x_fIII,y_fIII,z_fIII);
            fun2fIV=sph(nhs2,mhs2,x_fIV,y_fIV,z_fIV);
            fun2fV=sph(nhs2,mhs2,x_fV,y_fV,z_fV);
            fun2fVI=sph(nhs2,mhs2,x_fVI,y_fVI,z_fVI);
        elseif strcmp(axes,'oy')==1
            fun2fI=sph(nhs2,mhs2,x_fI,-z_fI,y_fI);
            fun2fII=sph(nhs2,mhs2,x_fII,-z_fII,y_fII);
            fun2fIII=sph(nhs2,mhs2,x_fIII,-z_fIII,y_fIII);
            fun2fIV=sph(nhs2,mhs2,x_fIV,-z_fIV,y_fIV);
            fun2fV=sph(nhs2,mhs2,x_fV,-z_fV,y_fV);
            fun2fVI=sph(nhs2,mhs2,x_fVI,-z_fVI,y_fVI);
        elseif strcmp(axes,'ox')==1
            fun2fI=sph(nhs2,mhs2,-z_fI,y_fI,x_fI);
            fun2fII=sph(nhs2,mhs2,-z_fII,y_fII,x_fII);
            fun2fIII=sph(nhs2,mhs2,-z_fIII,y_fIII,x_fIII);
            fun2fIV=sph(nhs2,mhs2,-z_fIV,y_fIV,x_fIV);
            fun2fV=sph(nhs2,mhs2,-z_fV,y_fV,x_fV);
            fun2fVI=sph(nhs2,mhs2,-z_fVI,y_fVI,x_fVI);
        end



        % APPROXIMATE SCALAR PRODUCT ON THE SPHERE
        [scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
            intsca_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
            fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI);

        NHS(i,j)=nhs2;
        MHS(i,j)=mhs2;
        scal(i,j)=scaf1f2;
    end
end

figure(1)
contourf(NHS,MHS,abs(scal));
xlabel('nhs')
ylabel('mhs')
colorbar

figure(2)
surf(NHS,MHS,abs(scal));
xlabel('nhs')
ylabel('mhs')
colorbar

figure(3)
surf(NHS,MHS,log10(abs(scal)));
xlabel('nhs')
ylabel('mhs')
title('log')
colorbar

figure(4)
contourf(NHS,MHS,log10(abs(scal).*(abs(scal)<10^-3)));
xlabel('nhs')
ylabel('mhs')
title('log')
grid minor
colorbar

figure(5)
plot_cs11(n,nn,real(fun1fI),real(fun1fII),real(fun1fIII),real(fun1fIV),real(fun1fV),real(fun1fVI))
title('ref. hs')

figure(6)
plot_cs11(n,nn,real(fun2fI),real(fun2fII),real(fun2fIII),real(fun2fIV),real(fun2fV),real(fun2fVI))
title('last hs tested')

fig_placier