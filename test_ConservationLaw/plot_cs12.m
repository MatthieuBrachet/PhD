function []=plot_cs12(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI,err)
% PLOT THE CUBED SPHERE GRID
% J-P. CROISILLE - JULY 20 2010
global radius

%% tracé de la courbe
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

surf(x1,y1,z1,funfI);hold on;
axis([-r r -r r -r r]);

% FACE E - II;
y2=sqrt(radius^2./delta);
x2=-xx.*y2;
z2=yy.*y2;

surf(x2,y2,z2,funfII);hold on;
axis([-r r -r r -r r]);

% FACE B - III;
x3=-sqrt(radius^2./delta);
y3=xx.*x3;
z3=-yy.*x3;

surf(x3,y3,z3,funfIII);hold on;
axis([-r r -r r -r r]);

% FACE W - IV;
y4=-sqrt(radius^2./delta);
x4=-xx.*y4;
z4=-yy.*y4;

surf(x4,y4,z4,funfIV);hold on;
axis([-r r -r r -r r]);

% FACE N - V;
z5=sqrt(radius^2./delta);
y5=xx.*z5;
x5=-yy.*z5;

surf(x5,y5,z5,funfV);hold on;
axis([-r r -r r -r r]);

% FACE S - VI;
z6=-sqrt(radius^2./delta);
y6=-xx.*z6;
x6=-yy.*z6;

surf(x6,y6,z6,funfVI);hold on;
axis([-r r -r r -r r]);

axis equal;
axis([-radius radius -radius radius -radius radius]);
shading interp
caxis([-err, err]);
colorbar;;

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

plot3(z,y,x,'k-','linewidth',2); hold on
plot3(z,-y,x,'k-','linewidth',2); hold on
plot3(z,y,-x,'k-','linewidth',2); hold on
plot3(z,-y,-x,'k-','linewidth',2); hold on

plot3(x,z,y,'k-','linewidth',2); hold on
plot3(-x,z,y,'k-','linewidth',2); hold on
plot3(x,z,-y,'k-','linewidth',2); hold on
plot3(-x,z,-y,'k-','linewidth',2); hold on

plotCoastLines(1,'k-');