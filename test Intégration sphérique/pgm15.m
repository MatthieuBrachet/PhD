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

N=64;
make_cs_grid(N);  

r=5;
nt=10;
nl=10;
[A,err_i]=compute_A_sym_alouges(r,nt,nl);
% eps_w=solve_weights_alouges(r,nt,nl,A,err_i,0,1);
weights=dxi*deta*dga;
% weights=dxi*deta*(dga+eps_w);

err_r=[];
for r=0:200,
    r
    err=[];
for theta=linspace(-pi/2,pi/2,10),    
    for lambda=linspace(-pi,pi,10),

        X(1)=r*cos(theta)*cos(lambda);
        X(2)=r*cos(theta)*sin(lambda);
        X(3)=r*sin(theta);
        funfI=fun_alouges(X,x_fI,y_fI,z_fI);
        funfII=fun_alouges(X,x_fII,y_fII,z_fII);
        funfIII=fun_alouges(X,x_fIII,y_fIII,z_fIII);
        funfIV=fun_alouges(X,x_fIV,y_fIV,z_fIV);
        funfV=fun_alouges(X,x_fV,y_fV,z_fV);
        funfVI=fun_alouges(X,x_fVI,y_fVI,z_fVI);
        [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
           int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);

        err=[err;abs(nrmg-4*pi*card_sin(r))];

    end
end
    err_r=[err_r;max(err)];
end
max(err)