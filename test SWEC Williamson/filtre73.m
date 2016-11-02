function [ ftr ] = filtre( na , opt_ftr )
global alfa_ftr
if strcmp(opt_ftr,'redonnet0') == 1
    f0=1; 
    f1=0;
    f2=0;
    f3=0;
    f4=0;
    f5=0;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=speye(size(ftr));
    ftri=speye(size(ftre));
    
elseif strcmp(opt_ftr,'redonnet2') == 1
    f0=1/2; 
    f1=1/4;
    f2=0;
    f3=0;
    f4=0;
    f5=0;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    ftri=speye(size(ftre));

elseif strcmp(opt_ftr,'redonnet4') == 1
    f0=10/16; 
    f1=4/16;
    f2=-1/16;
    f3=0;
    f4=0;
    f5=0;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    ftri=speye(size(ftre));
    
elseif strcmp(opt_ftr,'redonnet6') == 1
    f0=44/64; 
    f1=15/64;
    f2=-6/64;
    f3=1/64;
    f4=0;
    f5=0;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    ftri=speye(size(ftre));
    
    
elseif strcmp(opt_ftr,'redonnet8') == 1
    f0=186/256; 
    f1=56/256;
    f2=-28/256;
    f3=8/256;
    f4=-1/256;
    f5=0;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    ftri=speye(size(ftre));

elseif strcmp(opt_ftr,'redonnet10') == 1
    f0=772/1024; 
    f1=210/1024;
    f2=-120/1024;
    f3=45/1024;
    f4=-10/1024;
    f5=1/1024;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    ftri=speye(size(ftre));
    
elseif strcmp(opt_ftr,'visbal2') == 1
    aa=[0.5+alfa_ftr 0.5+alfa_ftr 0 0 0 0];
    f0=aa(1); 
    f1=aa(2)/2;
    f2=aa(3)/2;
    f3=aa(4)/2;
    f4=aa(5)/2;
    f5=aa(6)/2;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    
    ftri=sparse(eye(size(ftre))+alfa_ftr*diag(ones(length(ftre)-1,1),1)+alfa_ftr*diag(ones(length(ftre)-1,1),-1));
    ftri(1,end)=alfa_ftr; ftri(end,1)=alfa_ftr;
    
elseif strcmp(opt_ftr,'visbal4') == 1
    aa=[5/8+3/4*alfa_ftr 0.5+alfa_ftr -1/8+alfa_ftr/4 0 0 0];
    f0=aa(1); 
    f1=aa(2)/2;
    f2=aa(3)/2;
    f3=aa(4)/2;
    f4=aa(5)/2;
    f5=aa(6)/2;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    
    ftri=sparse(eye(size(ftre))+alfa_ftr*diag(ones(length(ftre)-1,1),1)+alfa_ftr*diag(ones(length(ftre)-1,1),-1));
    ftri(1,end)=alfa_ftr; ftri(end,1)=alfa_ftr;
    
elseif strcmp(opt_ftr,'visbal6') == 1
    aa=[11/16+5/8*alfa_ftr 15/32+17/16*alfa_ftr -3/16+3/8*alfa_ftr 1/32-1/16*alfa_ftr 0 0];
    f0=aa(1); 
    f1=aa(2)/2;
    f2=aa(3)/2;
    f3=aa(4)/2;
    f4=aa(5)/2;
    f5=aa(6)/2;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    
    ftri=sparse(eye(size(ftre))+alfa_ftr*diag(ones(length(ftre)-1,1),1)+alfa_ftr*diag(ones(length(ftre)-1,1),-1));
    ftri(1,end)=alfa_ftr; ftri(end,1)=alfa_ftr;
    
elseif strcmp(opt_ftr,'visbal8') == 1
    aa=[93/128+70/128*alfa_ftr 7/16+18/16*alfa_ftr -7/32+14/32*alfa_ftr 1/16-1/8*alfa_ftr -1/128+1/64*alfa_ftr 0];
    f0=aa(1); 
    f1=aa(2)/2;
    f2=aa(3)/2;
    f3=aa(4)/2;
    f4=aa(5)/2;
    f5=aa(6)/2;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    
    ftri=sparse(eye(size(ftre))+alfa_ftr*diag(ones(length(ftre)-1,1),1)+alfa_ftr*diag(ones(length(ftre)-1,1),-1));
    ftri(1,end)=alfa_ftr; ftri(end,1)=alfa_ftr;
    
elseif strcmp(opt_ftr,'visbal10') == 1
    aa=[193/256+126/256*alfa_ftr 105/256+302/256*alfa_ftr -15/64+30/64*alfa_ftr 45/512-90/512*alfa_ftr -5/256+10/256*alfa_ftr 1/512-2/512*alfa_ftr];
    f0=aa(1); 
    f1=aa(2)/2;
    f2=aa(3)/2;
    f3=aa(4)/2;
    f4=aa(5)/2;
    f5=aa(6)/2;
    
    %%
    lig1=[0,1, zeros(1,na-2)];
    col1=[zeros(na-1,1);1];
    sh1=toeplitz(col1,lig1);
    sh1i=inv(sh1);
    sh12=sh1^2;
    sh1i2=sh1i^2;
    sh13=sh12*sh1;
    sh1i3=sh1i2*sh1i;
    sh14=sh13*sh1;
    sh1i4=sh1i3*sh1i;
    sh15=sh14*sh1;
    sh1i5=sh1i4*sh1i;
    ftr=f0*eye(na)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
    ftre=sparse(ftr);
    
    ftri=sparse(eye(size(ftre))+alfa_ftr*diag(ones(length(ftre)-1,1),1)+alfa_ftr*diag(ones(length(ftre)-1,1),-1));
    ftri(1,end)=alfa_ftr; ftri(end,1)=alfa_ftr;
%    
end


ftr=ftri\ftre;

end

