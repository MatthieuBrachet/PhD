clc; clear all; close all;

N_cs=[1 2 3 4 5 8 10 12 14 18 20 25 30];
cs=6*N_cs.^2+2;
nhs_max=[3 6 9 12 15 24 30 34 40 52 58 73 86];
nb_hs=nhs_max.*(nhs_max+2);

figure(1)
plot(N_cs,nhs_max,'b-','Linewidth',2)
grid on
xlabel('N_cs')
ylabel('nhs_{max}')

figure(2)
plot(cs,nhs_max,'r-','Linewidth',2)
grid on
xlabel('Nb pts CS')
ylabel('nhs_{max}')

figure(3)
plot(cs,nb_hs,'g-','Linewidth',2)
grid on
xlabel('nb de pts sur la CS')
ylabel('nb de HS utiles (borne sup, il y en a trop)')

fig_placier

