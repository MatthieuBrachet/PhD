%%% COMPUTATION OF THE WEIGHTS (same for the 6 faces) FOR THE CUBED-SPHERE
%%% WITH THE PSEUDO-INVERSE METHOD OF FORNBERG-MARTEL
%%% => WITH EXACT PROPERTIES OF SYMMETRY (1/8 of the patch)

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

N=30;
make_cs_grid(N);  

[A,err_i]=compute_A_sym(0);

y=[];
kmin=1;
step=1;
kmax=size(A,1);
for k=kmin:step:kmax,
    
    k
    eps_w=solve_weights(A,err_i,k,1,1);
    y=[y;int_funs_fornberg(dxi*deta*(dga+eps_w))];
    
end

subplot(411);
p=plot(kmin:step:kmax,y(:,1));
p.LineWidth=1.5;
p.Color='black';
title('fonction constante =1');
subplot(412);
p=plot(kmin:step:kmax,y(:,2));
p.LineWidth=1.5;
p.Color='green';
title('f1');
subplot(413);
p=plot(kmin:step:kmax,y(:,3));
p.LineWidth=1.5;
p.Color='red';
title('f2');
subplot(414);
p=plot(kmin:step:kmax,y(:,4));
p.LineWidth=1.5;
p.Color='blue';
title('f3');
figure(1);
xlabel('k = nombre de HS prises en compte');
ylabel('erreur intégration');

x=[];
for i=1:length(y),
    x=[x;norm(y(i,:),2)];
end
[a,b]=min(x);
b
    
