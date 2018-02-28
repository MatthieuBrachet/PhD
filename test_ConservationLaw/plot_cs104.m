function []=plot_cs104(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
% plot solution on the cubed associated with the cs grid.
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;


%% plot solution
% Face IV
xi_IV=-x_fIV./y_fIV;
eta_IV=-z_fIV./y_fIV;
surf(xi_IV-2,eta_IV,funfIVe);
hold on
% Face I
xi_I=y_fI./x_fI;
eta_I=z_fI./x_fI;
surf(xi_I,eta_I,funfIe);
hold on
% Face II
xi_II=-x_fII./y_fII;
eta_II=z_fII./y_fII;
surf(xi_II+2,eta_II,funfIIe);
hold on
% Face III
xi_III=y_fIII./x_fIII;
eta_III=-z_fIII./x_fIII;
surf(xi_III+4,eta_III,funfIIIe);
hold on
% Face V
xi_V=-x_fV./z_fV;
eta_V=y_fV./z_fV;
surf(xi_V,eta_V+2,funfVe');
hold on
% Face VI
xi_VI=-x_fV./z_fV;
eta_VI=y_fV./z_fV;
surf(xi_VI,eta_VI-2,funfVIe');
hold on

%% CS edges
plot3(xi_I(1,:),eta_I(1,:),funfIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_I(end,:),eta_I(end,:),funfIe(end,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_I(:,end),eta_I(:,end),funfIe(:,end)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_I(:,1),eta_I(:,1),funfIe(:,1)+eps,'k','LineWidth',1.25); hold on;

plot3(xi_II(:,end)+2,eta_II(:,end),funfIIe(:,end)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_II(:,1)+2,eta_II(:,1),funfIIe(:,1)+eps,'k','LineWidth',1.25); hold on;

plot3(xi_III(1,:)+4,eta_III(1,:),funfIIIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_III(end,:)+4,eta_III(end,:),funfIIIe(end,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_III(:,end)+4,eta_III(:,end),funfIIIe(:,end)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_III(:,1)+4,eta_III(:,1),funfIIIe(:,1)+eps,'k','LineWidth',1.25); hold on;

plot3(xi_IV(:,end)-2,eta_IV(:,end),funfIVe(:,end)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_IV(:,1)-2,eta_IV(:,1),funfIVe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_IV(1,:)-2,eta_IV(1,:),funfIVe(1,:)+eps,'k','LineWidth',1.25); hold on;

plot3(xi_V(end,:),eta_V(end,:)+2,funfVe(:,end)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_V(:,1),eta_V(:,1)+2,funfVe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_V(:,end),eta_V(:,end)+2,funfVe(end,:)+eps,'k','LineWidth',1.25); hold on;

plot3(xi_VI(1,:),eta_VI(1,:)-2,funfVIe(:,1)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_VI(:,1),eta_VI(:,1)-2,funfVIe(1,:)+eps,'k','LineWidth',1.25); hold on;
plot3(xi_VI(:,end),eta_VI(:,end)-2,funfVIe(end,:)+eps,'k','LineWidth',1.25); hold on;



shading interp;
view(2)