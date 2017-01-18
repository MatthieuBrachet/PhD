%%% Courbes des maximum des |epsilon| en
%%% fonction de N

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
N_max=32;

y=[];

for N=4:2:N_max,
    
    N
    k=N^2/4;
    make_cs_grid(N);
    
    [A,err_i]=compute_A_sym(nhs_max);
    eps_w=solve_weights(A,err_i,k,1,1);
    y=[y;max(max(abs(eps_w)))];
    
end