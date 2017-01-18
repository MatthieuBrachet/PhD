%%% Test of the function from Alouges

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

N=8;
theta=0;
lambda=0;
X=zeros(3,1);
make_cs_grid(N); 
%k=N^2/4;
k=20;
[A,err_i]=compute_A_sym(0);
eps_w=solve_weights(A,err_i,k,1,1);
% eps_w=zeros(nn,nn);
weights=dxi*deta*(dga+eps_w);

y=[];
r=20;
for r=1:1:200,
    err=[];
    r
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

    err=[err;abs(nrmg-4*pi*sin(norm(X,2))/norm(X,2))];
    
end
end
    y=[y;max(err)];
end

% x=[x_fI,x_fII,x_fIII,x_fIV,x_fV,x_fVI];
% y=[y_fI,y_fII,y_fIII,y_fIV,y_fV,y_fVI];
% z=[z_fI,z_fII,z_fIII,z_fIV,z_fV,z_fVI];
% fun=[funfI,funfII,funfIII,funfIV,funfV,funfVI];
% mesh(x,y,z,real(fun));
% colorbar;
% r=1;
% axis([-r r -r r -r r]);
% axis equal;
