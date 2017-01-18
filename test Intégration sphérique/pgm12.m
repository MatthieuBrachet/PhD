%%% Courbes d'erreurs des formules avec et sans les epsilon
%%% en fonction de N, pour les H.S.

clear all;
global n nn;
global radius;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=4;
make_cs_grid(N);
nhs_max=0;
nhs=40;
N_max=32;

y=[];
z=[];

for N=4:2:N_max,
    
    N
    k=N^2/4;
    make_cs_grid(N);
    
    weights=dxi*deta*dga;
    tmp=view_int_sph(nhs,weights,1,1);
    y=[y;tmp(4)];
    
    [A,err_i]=compute_A_sym(nhs_max);
    eps_w=solve_weights(A,err_i,k,1,1);
    weights=dxi*deta*(dga+eps_w);
    tmp=view_int_sph(nhs,weights,1,1);
    z=[z;tmp(4)];
    
end

p=plot(log(4:2:N_max),log(y));
p.Marker='x';
hold on;
pp=plot(log(4:2:N_max),log(z));
pp.Marker='*';
pp.LineWidth=1.5;
grid minor; grid on; title('Erreurs en fonction de N, pour (n,m)=(40,16)');
xlabel('log(N)'); ylabel('log(erreur)');




