clear all;
global n nn;
global mm na nb;
global radius;
global xi eta dxi deta xx yy delta deltab;
global alfa beta;
global alfacr betacr;
global alfa1;
global alfag betag;
global p k;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global p1 k1;
global  aaa bbb itestop
%%%%%%%%%%%%%%%%%%%
global alphad u0; % CF fun4.m, option 9= cas test 1 de Williamson;
%
mod99; % APPEL MODULE "PROBLEME"
%
% ASSEMBLING THE VECTOR FUNCTION TO TAKE THE DIVERGENCE
% ON FACES I, II, III, IV.
% 1-  DIVERGENCE OF A VECTOR FIELDS GIVEN IN FUN6.
% c_fI=[0 0 0];
fun_fI=ones(n+1,n+1);
fun_fII=2*ones(n+1,n+1);
fun_fIII=3*ones(n+1,n+1);
fun_fIV=4*ones(n+1,n+1);
fun_fV=5*ones(n+1,n+1);
fun_fVI=6*ones(n+1,n+1);
%
fig=figure(10);clf(fig);
% pos=get(fig,'Position');
% pos([3 4])=400;
% set(fig,'Position',pos);
% FigHandle = figure;
set(fig, 'Position', [100, 100, 1049, 895]);
% ax=axes('Parent', fig, 'pos',[0 0 1 1]);
ax=axes('Parent', fig);
colormap('gray');
set(gcf,'color','white');
plot_cs15(n,nn,fun_fI,fun_fII,fun_fIII,fun_fIV,fun_fV,fun_fVI);colorbar;
colorbar('off');
axis(ax,'equal');
axis(ax,'off');
% view(ax,3);
view(ax,135,35);
% title('CUBED SPHERE'); 
hold on;
%%% WRITE TEXT. VALUE OF A1(1) IN THE NEXT LINE= ESTIMATED CONVERGENCE RATE FOR
%%% PSI. ALSO ADJUST POSITION OF THE TEXT
% annotation('Panel (I)',[0 100], [0 100]);
% th=annotation('textarrow',[.3,.6],[.7,.4],'String','ABC');
th=annotation('textarrow',[.15,.4],[.1,.35],'String','PANEL (I)','HeadStyle','cback2');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
th=annotation('textarrow',[.85,.65],[.1,.36],'String','PANEL (II)');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
th=annotation('textarrow',[.85,.71],[.85,.64],'String','PANEL (III)');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
th=annotation('textarrow',[.15,0.32],[.85,.65],'String','PANEL (IV)');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
th=annotation('textarrow',[.65,0.60],[.95,.72],'String','PANEL (V)');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
th=annotation('textarrow',[.5,0.51],[.05,.26],'String','PANEL (VI)');
set(th,'FontSize',14);
set(th,'FontWeight','demi');
% figure(11);
% xp1=1;yp1=0;zp1=0;
% vx1=1;vy1=0;vz1=0;
% annotation('arrow',xp,yp,zp);
% scale=1;
% hq=quiver3(xp1,yp1,zp1,vx1,vy1,vz1,scale);
% set(hq,'LineWidth',2,'Color','k','MaxHeadSize',4);
% thx=annotation('textbox',[.15,.3,0.1,0.1],'String','X');
% thx=annotation('text',[0.1,0.1],'String','X');
xp1=[1,2];
yp1=[0,0];
zp1=[0,0];
hlx=line(xp1,yp1,zp1,'LineWidth',3,'Color','k');
text(2.1,0.,0.,'X','FontSize',16);
%%%%%%
xp2=[0,0];
yp2=[1,2];
zp2=[0,0];
hly=line(xp2,yp2,zp2,'LineWidth',3,'Color','k');
text(0.,2.1,0.,'Y','FontSize',16);
%%%%%%
xp3=[0,0];
yp3=[0,0];
zp3=[1,2];
hly=line(xp3,yp3,zp3,'LineWidth',3,'Color','k');
text(0.,0.,2.1,'Z','FontSize',16);
%%%
% text(50,50,'X','FontSize',14);
% text(10,10,'X');
% set(gcf,'LineWidth',20); hold on;
%ht=text('Position',[-10 -10],'String','slope for \psi: 5.66');
% ht=text('Position',[-1.3,-6.5,0],'String','slope for \psi: 5.66');

