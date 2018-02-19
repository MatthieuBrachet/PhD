global n X Y
global opt_ftr
global gp hp omega coriolis
global FTR kx ky px py

omega=7.292d-05;
teta0=pi/2; coriolis=2.*omega.*sin(teta0);
gp=9.80616;
hp=10;

h=1./n;
x=h:h:1;
[X,Y]=meshgrid(x,x);
X=reshape(X,[],1);
Y=reshape(Y,[],1);

id=speye(n,n);
k=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
k(1,end)=-1; k(end,1)=1;
k=sparse(k./(2*h));
p=1/6*diag(ones(n-1,1),1)+1/6*diag(ones(n-1,1),-1)+4/6*eye(n,n);
p(end,1)=1/6; p(1,end)=1/6;
p=sparse(p);

kx=kron(id,k);
ky=kron(k,id);
px=kron(id,p);
py=kron(p,id);


if strcmp(opt_ftr,'redonnet10')==1
    ftr0=772/1024;
    ftr1=420/1024;
    ftr2=-240/1024;
    ftr3=90/1024;
    ftr4=-20/1024;
    ftr5=2/1024;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet8')==1
    ftr0=186/256; 
    ftr1=112/256;
    ftr2=-56/256;
    ftr3=16/256;
    ftr4=-2/256;
    ftr5=0;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet6')==1
    ftr0=44/64;
    ftr1=30/64;
    ftr2=-12/64;
    ftr3=2/64;
    ftr4=0;
    ftr5=0;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2))+ftr3/2*(J^3+J^(n-3))+ftr4/2*(J^4+J^(n-4))+ftr5/2*(J^5+J^(n-5));
    
elseif strcmp(opt_ftr,'redonnet4')==1
    ftr0=10/16;
    ftr1=8/16;
    ftr2=-2/16;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1))+ftr2/2*(J^2+J^(n-2));
elseif strcmp(opt_ftr,'redonnet2')==1
    ftr0=1/2;
    ftr1=1/2;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1/2*(J+J^(n-1));
elseif strcmp(opt_ftr,'bogey6')==1
    d0=0.234810479761700;
    d1=-.199250131285813;
    d2=0.120198310245186;
    d3=-.049303775636020;
    d4=0.012396449873964;
    d5=-.001446093078167;
    
    ftr0=1-d0;
    ftr1=-d1;
    ftr2=-d2;
    ftr3=-d3;
    ftr4=-d4;
    ftr5=-d5;
    J=diag(ones(n-1,1),1); J(end,1)=1;
    ftr=ftr0.*eye(n,n)+ftr1*(J+J^(n-1))+ftr2*(J^2+J^(n-2))+ftr3*(J^3+J^(n-3))+ftr4*(J^4+J^(n-4))+ftr5*(J^5+J^(n-5));
else
    ftr=speye(n,n);
end
ftr=sparse(ftr);

ftrx=kron(id,ftr);
ftry=kron(ftr,id);
FTR=ftrx*ftry;

