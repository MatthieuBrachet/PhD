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
%figure(1);
N=2;
n=N-1;nn=N+1;
%plot_cs1_mesh(n,nn);
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
x_cs1=zeros(kN,1);y_cs1=zeros(kN,1);z_cs1=zeros(kN,1);
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
x_cs1(1:8)=x_cs(1:8);
y_cs1(1:8)=y_cs(1:8);
z_cs1(1:8)=z_cs(1:8);
%
ikomp=8;
% b) 12 (N-1) edge points
% bottom edges face I
for k=1:N-1;
    x_cs1(ikomp+k)=x_fI(k+1,1);
    y_cs1(ikomp+k)=y_fI(k+1,1);
    z_cs1(ikomp+k)=z_fI(k+1,1);
end
ikomp=ikomp+N-1;
% bottom edges face II
for k=1:N-1;
    x_cs1(ikomp+k)=x_fII(k+1,1);
    y_cs1(ikomp+k)=y_fII(k+1,1);
    z_cs1(ikomp+k)=z_fII(k+1,1);
end
ikomp=ikomp+N-1;
% bottom edges face III
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIII(k+1,1);
    y_cs1(ikomp+k)=y_fIII(k+1,1);
    z_cs1(ikomp+k)=z_fIII(k+1,1);
end
ikomp=ikomp+N-1;
% bottom edges face IV
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIV(k+1,1);
    y_cs1(ikomp+k)=y_fIV(k+1,1);
    z_cs1(ikomp+k)=z_fIV(k+1,1);
end
ikomp=ikomp+N-1;
%%%%%%%%%%%%%%%%%%%%%%%%%
% top edges face I
for k=1:N-1;
    x_cs1(ikomp+k)=x_fI(k+1,nn);
    y_cs1(ikomp+k)=y_fI(k+1,nn);
    z_cs1(ikomp+k)=z_fI(k+1,nn);
end
ikomp=ikomp+N-1;
%%%%%%%%%%%%%%%%%%%%%%%%%
% top edges face II
for k=1:N-1;
    x_cs1(ikomp+k)=x_fII(k+1,nn);
    y_cs1(ikomp+k)=y_fII(k+1,nn);
    z_cs1(ikomp+k)=z_fII(k+1,nn);
end
ikomp=ikomp+N-1;
%%%%%%%%%%%%%%%%%%%%%%%%
% top edges face III
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIII(k+1,nn);
    y_cs1(ikomp+k)=y_fIII(k+1,nn);
    z_cs1(ikomp+k)=z_fIII(k+1,nn);
end
ikomp=ikomp+N-1;
%%%%%%%%%%%%%%%%%%%%%%%%
% top edges face IV
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIV(k+1,nn);
    y_cs1(ikomp+k)=y_fIV(k+1,nn);
    z_cs1(ikomp+k)=z_fIV(k+1,nn);
end
ikomp=ikomp+N-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% side edges face I
for k=1:N-1;
    x_cs1(ikomp+k)=x_fI(1,k+1);
    y_cs1(ikomp+k)=y_fI(1,k+1);
    z_cs1(ikomp+k)=z_fI(1,k+1);
end
ikomp=ikomp+N-1;
% side edges face II
for k=1:N-1;
    x_cs1(ikomp+k)=x_fII(1,k+1);
    y_cs1(ikomp+k)=y_fII(1,k+1);
    z_cs1(ikomp+k)=z_fII(1,k+1);
end
ikomp=ikomp+N-1;
% side edges face III
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIII(1,k+1);
    y_cs1(ikomp+k)=y_fIII(1,k+1);
    z_cs1(ikomp+k)=z_fIII(1,k+1);
end
ikomp=ikomp+N-1;
% side edges face IV
for k=1:N-1;
    x_cs1(ikomp+k)=x_fIV(1,k+1);
    y_cs1(ikomp+k)=y_fIV(1,k+1);
    z_cs1(ikomp+k)=z_fIV(1,k+1);
end
ikomp=ikomp+N-1;
%
% ikomp
% fprintf('-----\n');
%%%%%%%%%%%%%%%%%%%%%%%%%
% c) interior points of panels
% interior points of face I
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fI(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fI(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fI(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
% interior points of face II
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fII(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fII(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fII(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
% interior points of face III
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fIII(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fIII(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fIII(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
% interior points of face IV
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fIV(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fIV(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fIV(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
% interior points of face V
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fV(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fV(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fV(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
% interior points of face VI
for k2=1:N-1;
    for k1=1:N-1;
        k=(k2-1)*(N-1)+k1;
        x_cs1(ikomp+k)=x_fVI(1+k1,1+k2);
        y_cs1(ikomp+k)=y_fVI(1+k1,1+k2);
        z_cs1(ikomp+k)=z_fVI(1+k1,1+k2);
    end
end
ikomp=ikomp+(N-1)^2;
co_cs1=[x_cs1,y_cs1,z_cs1];
co_cs=co_cs1;
x_cs=x_cs1;y_cs=y_cs1;z_cs=z_cs1;
ikomp
%break;
% aa=zeros(9,26);
%
nhs_max=3;
nhs_maxd=(nhs_max+1)^2;
%aat=zeros(nhs_maxd,26);
aab=[];
for ikomp=1:nhs_max+1,
    nhs2=ikomp-1;
    aan=zeros(2*nhs2+1,26);
    mhs2=0;
    fun_cs00=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
    for pl=1:kN,
       aan(1,pl)=fun_cs00(pl);
    end
    for mhs2=1:nhs2,
      fun_csa=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
      for pl=1:kN,
         aan(1+(2*mhs2-1),pl)=conj(fun_csa(pl));
      end
      for pl=1:kN,
        aan(2+(2*mhs2-1),pl)=fun_csa(pl);
      end          
    end
    rank(aan)
    aab=[aab;aan];
    rank(aab)
end
%break;
fprintf('-------------\n');
aa0=zeros(1,26);
nhs2=0;mhs2=0;
fun_cs00=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa0(1,pl)=fun_cs00(pl);
end
%
rank(aa0)
aa=aa0;
rank(aa)
%
aa1=zeros(3,26);
nhs2=1;mhs2=0;
fun_cs10=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa1(1,pl)=fun_cs10(pl);
end
nhs2=1;mhs2=1;
fun_cs11=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa1(2,pl)=conj(fun_cs11(pl));
end
for pl=1:kN,
    aa1(3,pl)=fun_cs11(pl);
end
%%%%%%%%%
rank(aa1)
size(aa1)
aa=[aa;aa1];
rank(aa)
size(aa)
% aa1=aa(1:4,1:26);
% rank(aa1)
%%%%%%%%%
aa2=zeros(5,26);
nhs2=2;mhs2=0;
fun_cs20=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa2(1,pl)=fun_cs20(pl);
end
nhs2=2;mhs2=1;
fun_cs21=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa2(2,pl)=conj(fun_cs21(pl));
end
for pl=1:kN,
    aa2(3,pl)=fun_cs21(pl);
end
nhs2=2;mhs2=2;
fun_cs22=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa2(4,pl)=conj(fun_cs22(pl));
end
for pl=1:kN,
    aa2(5,pl)=fun_cs22(pl);
end
%%%%%%%%%%
rank(aa2)
aa=[aa;aa2];
rank(aa)
% aa2=aa(1:9,1:26);
% rank(aa2)
%%%%%%%%%%
aa3=zeros(7,26);
nhs2=3;mhs2=0;
fun_cs30=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa3(1,pl)=fun_cs30(pl);
end
nhs2=3;mhs2=1;
fun_cs31=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa3(2,pl)=conj(fun_cs31(pl));
end
for pl=1:kN,
    aa3(3,pl)=fun_cs31(pl);
end
nhs2=3;mhs2=2;
fun_cs32=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa3(4,pl)=conj(fun_cs32(pl));
end
for pl=1:kN,
    aa3(5,pl)=fun_cs32(pl);
end
nhs2=3;mhs2=3;
fun_cs33=sph(nhs2,mhs2,x_cs,y_cs,z_cs);
for pl=1:kN,
    aa3(6,pl)=conj(fun_cs33(pl));
end
for pl=1:kN,
    aa3(7,pl)=fun_cs33(pl);
end
%%%%%%%%%%
%%%%%%%%%%
rank(aa3)
aa=[aa;aa3];
rank(aa)
%%%%%%%%%%
% aa3=zeros(7,26);
% aa3=aa(10:16,1:26);
% rank(aa3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         disp(['n''       = ' num2str(nhs2)])
%         disp(['m''       = ' num2str(mhs2)])
%         disp(['scal     = ' num2str(scaf1f2)])
%         disp(['abs(scal) = ' num2str(abs(scaf1f2))])
%         disp('**********************')
%     end
% end