clc; clear all; close all;
format long;

N=25;
n=N-1;nn=N+1;
make_cs_grid(N);
eps=0.5;

nhs_max=6*N^2+2;

harm=[];
dim=[];
iter=[];
rg=-1;
nhs=0;
while nhs<=nhs_max && rg<nhs_max
    if rg<nhs_max
        mhs=-nhs;
        while mhs<=nhs && rg<=nhs_max
            [ hs ] = sph_cs( nhs, mhs );
            harm=[harm hs];
            mhs=mhs+1;
        end
        rg=rank(harm,eps);
    else 
        rg=nhs_max;
    end
    iter=[iter nhs];
    dim=[dim rg];
    clc; disp(num2str(rg/nhs_max*100))
    nhs=nhs+1;
end
I=ones(size(dim));

figure(1)
plot(iter,dim,'x-',iter,nhs_max.*I,'Linewidth',2)
grid on
legend('dim(HS) sur la CS','dim(CS)')
xlabel('nhs')
ylabel('rang')

figure(2)
plot(iter,[1 diff(dim)],'x-','Linewidth',2)
grid on
ylabel('gain sur le rang')
xlabel('nhs')

figure(3)
plot(iter,(dim)./((iter+1).^2),'x-','Linewidth',2)
grid on
xlabel('nhs')
ylabel('dim/nhs^2')

fig_placier

iter(end)