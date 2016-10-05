function []=plot_cs16(n)
% PLOT THE CUBED SPHERE GRID
% present the C-S grid with coasts.
global radius
radius = 1;
nn=n+2;
%% tracé du fond
r=radius;
%
xi=linspace(-pi/4, pi/4, nn); 
eta=linspace(-pi/4, pi/4, nn); 
for i=1:nn,
  for j=1:nn,
    xx(i,j)=tan(xi(i));
    yy(i,j)=tan(eta(j));
  end
end
%
for i=1:nn,
  for j=1:nn,
    delta(i,j)=1+xx(i,j)^2+yy(i,j)^2;
  end
end

% FACE F - I;
x1=sqrt(radius^2./delta);
y1=xx.*x1;
z1=yy.*x1;

fun=ones(size(x1));
colormap([0.9 0.9 0.9]); hold on
axis equal
set(gca,'Visible','off');
surf(x1,y1,z1,fun);hold on;

% FACE E - II;
y2=sqrt(radius^2./delta);
x2=-xx.*y2;
z2=yy.*y2;

surf(x2,y2,z2,fun);hold on;

% FACE B - III;
x3=-sqrt(radius^2./delta);
y3=xx.*x3;
z3=-yy.*x3;

surf(x3,y3,z3,fun);hold on;

% FACE W - IV;
y4=-sqrt(radius^2./delta);
x4=-xx.*y4;
z4=-yy.*y4;

surf(x4,y4,z4,fun);hold on;

% FACE N - V;
z5=sqrt(radius^2./delta);
y5=xx.*z5;
x5=-yy.*z5;

surf(x5,y5,z5,fun);hold on;

% FACE S - VI;
z6=-sqrt(radius^2./delta);
y6=-xx.*z6;
x6=-yy.*z6;

surf(x6,y6,z6,fun);hold on;

%% add the borders of each panel
[~, ts,~]=cart2sph(x1(1), y1(1), z1(1));
ts=abs(ts);

th=linspace(-ts,ts,n+1);
pp=pi/4*ones(1,n+1);

[x,y,z]=sph2cart(pp,th,radius);
plot3(x,y,z,'k-','linewidth',2); hold on
plot3(-x,y,z,'k-','linewidth',2); hold on
plot3(x,-y,z,'k-','linewidth',2); hold on
plot3(-x,-y,z,'k-','linewidth',2); hold on
plot3(x1,y1,z1,'k-','linewidth',0.5);

plot3(z,y,x,'k-','linewidth',2); hold on
plot3(z,-y,x,'k-','linewidth',2); hold on
plot3(z,y,-x,'k-','linewidth',2); hold on
plot3(z,-y,-x,'k-','linewidth',2); hold on

plot3(x,z,y,'k-','linewidth',2); hold on
plot3(-x,z,y,'k-','linewidth',2); hold on
plot3(x,z,-y,'k-','linewidth',2); hold on
plot3(-x,z,-y,'k-','linewidth',2); hold on


%% add the mesh
plot3(x1,y1,z1,'k-','linewidth',0.5);hold on;
plot3(x1',y1',z1','k-','linewidth',0.5);hold on;
plot3(x2,y2,z2,'k-','linewidth',0.5);hold on;
plot3(x2',y2',z2','k-','linewidth',0.5);hold on;
plot3(x3,y3,z3,'k-','linewidth',0.5);hold on;
plot3(x3',y3',z3','k-','linewidth',0.5);hold on;
plot3(x4,y4,z4,'k-','linewidth',0.5);hold on;
plot3(x4',y4',z4','k-','linewidth',0.5);hold on;
plot3(x5,y5,z5,'k-','linewidth',0.5);hold on;
plot3(x5',y5',z5','k-','linewidth',0.5);hold on;
plot3(x6,y6,z6,'k-','linewidth',0.5);hold on;
plot3(x6',y6',z6','k-','linewidth',0.5);hold on;


%% add the coast
plotCoastLines(1,'b-');