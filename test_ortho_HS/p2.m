clc; clear all; close all;
format long;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
N=1;
n=N-1;nn=N+1;
plot_cs1_mesh(n,nn);
make_cs_grid(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear system for the space Y_n = HS
% of degree n
% nhs1=1;
% nhsz=2*nhs1+1;
kN=6*N^2+2;
%
%aa=zeros(nhsz,kN);
%
x_cs=zeros(kN,1);y_cs=zeros(kN,1);z_cs=zeros(kN,1);
% listing the ponts of the CS
% 1) the 8 vertices of the CS
x_cs(1)=x_fI(1,1);
y_cs(1)=y_fI(1,1);
z_cs(1)=z_fI(1,1);
%
x_cs(2)=x_fI(nn,1);
y_cs(2)=y_fI(nn,1);
z_cs(2)=z_fI(nn,1);
%
x_cs(3)=x_fI(1,nn);
y_cs(3)=y_fI(1,nn);
z_cs(3)=z_fI(1,nn);
%
x_cs(4)=x_fI(nn,nn);
y_cs(4)=y_fI(nn,nn);
z_cs(4)=z_fI(nn,nn);
%
x_cs(5)=x_fIII(1,1);
y_cs(5)=y_fIII(1,1);
z_cs(5)=z_fIII(1,1);
%
x_cs(6)=x_fIII(nn,1);
y_cs(6)=y_fIII(nn,1);
z_cs(6)=z_fIII(nn,1);
%
x_cs(7)=x_fIII(1,nn);
y_cs(7)=y_fIII(1,nn);
z_cs(7)=z_fIII(1,nn);
%
x_cs(8)=x_fIII(nn,nn);
y_cs(8)=y_fIII(nn,nn);
z_cs(8)=z_fIII(nn,nn);
%
aa=zeros(9,8);
%
nhs2=0;mhs2=0;
fun_cs00=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(1,pl)=fun_cs00(pl);
end
nhs2=1;mhs2=0;
fun_cs10=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(2,pl)=fun_cs10(pl);
end
nhs2=1;mhs2=1;
fun_cs11=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(3,pl)=conj(fun_cs11(pl));
end
for pl=1:kN,
    aa(4,pl)=fun_cs11(pl);
end
rank(aa)
nhs2=2;mhs2=0;
fun_cs20=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(5,pl)=fun_cs20(pl);
end
nhs2=2;mhs2=1;
fun_cs21=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(6,pl)=conj(fun_cs21(pl));
end
for pl=1:kN,
    aa(7,pl)=fun_cs21(pl);
end
nhs2=2;mhs2=2;
fun_cs22=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa(8,pl)=conj(fun_cs22(pl));
end
for pl=1:kN,
    aa(9,pl)=fun_cs22(pl);
end
rank(aa)
aa1=zeros(5,8);
aa1=aa(5:9,1:8);
rank(aa1)

%         
%         disp(['n''       = ' num2str(nhs2)])
%         disp(['m''       = ' num2str(mhs2)])
%         disp(['scal     = ' num2str(scaf1f2)])
%         disp(['abs(scal) = ' num2str(abs(scaf1f2))])
%         disp('**********************')
%     end
% end