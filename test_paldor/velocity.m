function [ C ] = velocity(hovmoller,ddt)
[Ntps,Nspace]=size(hovmoller);
for ii = 1 : Ntps-1
    iter=[1:Nspace-1];
    ff1=fft(hovmoller(ii,:));
    ff2=fft(hovmoller(ii+1,:));
    angular=angle(ff2(2:end)./ff1(2:end));
    velocity=-Nspace./(2*pi*iter*ddt).*angular;
    vel(ii)=mean(velocity);
end
C=1./mean(vel);