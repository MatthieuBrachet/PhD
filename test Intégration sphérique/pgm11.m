%%% Method in the book of Freeden et al. p.1201, part 2.2.
%%% to have a measure of the (geometric) quality of the chosen grid

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
v=[];
v1=[];
v2=[];
n=10000;
phi=linspace(-pi/4,pi/4,n);
theta=linspace(-pi/4,pi/4,n);
[x,y,z]=sph2cart(theta,phi,1);
space=[];
for i=1:n,
    p=[x(i),y(i),z(i)];
    if p(1)>=x_fI(1,1) && p(1)<=x_fI(N/2+1,N/2+1) && ...
       p(2)>=y_fI(1,N/2+1) && p(2)<=y_fI(N+1,N/2+1) && ...
       p(3)>=z_fI(N/2+1,1) && p(3)<=z_fI(N/2+1,N+1),
        space=[space;p];
    end
end

for N=4:2:32,
    
    N
    make_cs_grid(N);
    points=[reshape(x_fI,[(N+1)^2,1]),reshape(y_fI,[(N+1)^2,1]),reshape(z_fI,[(N+1)^2,1])];
    dist1=[];
    for i=1:size(space,1),
        dist1=[dist1;dist_to_set(space(i,:),points)];
    end

    dist2=[];
    for i=1:size(points,1),
        dist2=[dist2;dist_to_set(points(i,:),[points(1:i-1,:);points(i+1:end,:)])];
    end

    v1=[v1;max(dist1)];
    v2=[v2;min(dist2)];
    v=[v;2*max(dist1)/min(dist2)];
    
end
            
subplot(311);
p=plot(4:2:32,v1);
p.LineWidth=1.5;
p.Marker='x';
p.Color='black';
title('h(X_N) en fonction de N');
ylabel('h(X_N)');
subplot(312);
p=plot(4:2:32,v2);
p.LineWidth=1.5;
p.Marker='x';
p.Color='green';
title('delta(X_N) en fonction de N');
ylabel('delta(X_N)');
subplot(313);
p=plot(4:2:32,v);
p.LineWidth=1.5;
p.Marker='x';
p.Color='red';
title('ratio de grille en fonction de N');
ylabel('ratio');
xlabel('N');
            
            
            
            
            
            
            
