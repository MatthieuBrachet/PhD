function [ ftr ] = filtre( na , opt_ftr )

if opt_ftr == 0
    f0=1; 
    f1=0;
    f2=0;
    f3=0;
    f4=0;
    f5=0;
    
elseif opt_ftr == 2
    f0=1/2; 
    f1=1/4;
    f2=0;
    f3=0;
    f4=0;
    f5=0;

elseif opt_ftr == 4
    f0=10/16; 
    f1=4/16;
    f2=-1/16;
    f3=0;
    f4=0;
    f5=0;
    
elseif opt_ftr == 6
    f0=44/64; 
    f1=15/64;
    f2=-6/64;
    f3=1/64;
    f4=0;
    f5=0;
    
elseif opt_ftr == 8
    f0=186/256; 
    f1=56/256;
    f2=-28/256;
    f3=8/256;
    f4=-1/256;
    f5=0;

elseif opt_ftr == 10
    f0=772/1024; 
    f1=210/1024;
    f2=-120/1024;
    f3=45/1024;
    f4=-10/1024;
    f5=1/1024;
end



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
ftr=sparse(ftr);
end

