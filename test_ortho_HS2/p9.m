clc; clear all; close all;
format long;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

tol=.01;
N=5;
n=N-1;nn=N+1;
make_cs_grid(N);

kN=6*N^2+2;  
nhs_max=kN;

for nhs=0:nhs_max
    A=[];
    for mhs=-nhs:nhs
        clc; [nhs mhs]
        res_fI = sph( nhs,mhs,x_fI,y_fI,z_fI );
        vect=reshape(res_fI,1,[]);
        res_fII = sph( nhs,mhs,x_fII,y_fII,z_fII );
        vect=[vect reshape(res_fII,1,[])];
        res_fIII = sph( nhs,mhs,x_fIII,y_fIII,z_fIII );
        vect=[vect reshape(res_fIII,1,[])];
        res_fIV = sph( nhs,mhs,x_fIV,y_fIV,z_fIV );
        vect=[vect reshape(res_fIV,1,[])];
        res_fV = sph( nhs,mhs,x_fV,y_fV,z_fV );
        vect=[vect reshape(res_fV,1,[])];
        res_fVI = sph( nhs,mhs,x_fVI,y_fVI,z_fVI );
        vect=[vect reshape(res_fVI,1,[])];
        
        A=[A; vect];
    end
    rang(nhs+1)=rank(A*A',tol);
end

figure(1)
plot(rang/kN)