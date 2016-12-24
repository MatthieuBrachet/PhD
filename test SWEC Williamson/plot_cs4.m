function []=plot_cs4(n,nn)
% images de panels
% PLOT THE CUBED SPHERE GRID
% J-P. CROISILLE - JULY 20 2010
% clear all;
% n=3,nn=5; 
% n=7, nn=9; % number of points by face
% n=15, nn=17; % number of points by face
%  n=31, nn=33; % number of points by face
% n =63, nn=65; % number of points by face
 % n=127, nn=129; % number of points by face
% n=255, nn=257; % number of points by face
% n=511, nn=513; % number of points by face
%
radius=1;r=radius;
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
%
%% PLOT OF THE BOUNDARY OF THE ALFA/BETA DOMAIN FOR ONE FACE
%%
alfaL=-atan(cos(eta));
betaL=atan(tan(eta)/sqrt(2));
%
alfaR=atan(cos(eta));
betaR=atan(tan(eta)/sqrt(2));
%
alfaB=atan(tan(xi)/sqrt(2));
betaB=-atan(cos(xi));
%
alfaT=atan(tan(xi)/sqrt(2));
betaT=atan(cos(xi));
%
plot(alfaL,betaL,alfaR,betaR,alfaB,betaB,alfaT,betaT); grid;
%
alfa=zeros(nn,nn);
beta=zeros(nn,nn);
%
for i=1:nn,
  for j=1:nn,
        alfa(i,j)=atan(tan(xi(i))/sqrt(1+tan(eta(j))^2));
        beta(i,j)=atan(tan(eta(j))/sqrt(1+tan(xi(i))^2));
  end
end
%
% 2- PLOT OF A PATCH OF THE CUBED SPHERE. REPRESENTATION DES NAGLES ALFA/BETA
% LIEN AVEC LA REPRESENTATION "EN PERSPECTIVE" DE RONCCHI PAR ETABLIE
% CLAIREMENT.
figure(2);
% PLOT HORIZONTAL ISOLINES
for j=1:nn,
    plot(alfa(:,j),beta(:,j)); hold on;
end
% PLOT VERTICAL ISOLINES
for i=1:nn,
    plot(alfa(i,:),beta(i,:)); hold on;
end
axis([-1 1 -1 1]);
title('A PATCH OF THE CUBED SPHERE. N=15');
axis equal;
% 3- DESSIN DES 6 CARRES POUR ARTCILE CS1.
% un carre
figure(3);
xc=zeros(2,2);yc=zeros(2,2);
for j=1:2,
        xc(:,j)=0:1;
        yc(j,:)=0:1;
end
% FACE I
for i=1:2,
    plot(xc(i,:),yc(i,:),'-'); hold on;
    plot(xc(:,i),yc(:,i),'-'); hold on;
end
text(0.5,0.5,'F(I)');
% FACE II
for i=1:2,
    plot(xc(i,:)+1,yc(i,:),'-'); hold on;
    plot(xc(:,i)+1,yc(:,i),'-'); hold on;
end
text(0.5+1,0.5,'E(II)');
% FACE III
for i=1:2,
    plot(xc(i,:)+2,yc(i,:),'-'); hold on;
    plot(xc(:,i)+2,yc(:,i),'-'); hold on;
end
text(0.5+2,0.5,'B(III)');
% FACE IV
for i=1:2,
    plot(xc(i,:)-1,yc(i,:),'-'); hold on;
    plot(xc(:,i)-1,yc(:,i),'-'); hold on;
end
text(0.5-1,0.5,'W(IV)');
% FACE V
for i=1:2,
    plot(xc(i,:),yc(i,:)+1,'-'); hold on;
    plot(xc(:,i),yc(:,i)+1,'-'); hold on;
end
text(0.5,0.5+1,'N(V)');
% FACE VI
for i=1:2,
    plot(xc(i,:),yc(i,:)-1,'-'); hold on;
    plot(xc(:,i),yc(:,i)-1,'-'); hold on;
end
text(0.5,0.5-1,'S(VI)');
axis([-1.5 3.5 -1.5 2.5]);
title('NOTATION FOR THE 6 PATCHES OF THE CUBED-SPHERE');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3- DESSIN DES 6 CARRES POUR ARTCILE CS1.
% un carre
figure(4);
xc=zeros(2,2);yc=zeros(2,2);
for j=1:2,
        xc(:,j)=0:1;
        yc(j,:)=0:1;
end
% FACE I
for i=1:2,
    plot(xc(i,:),yc(i,:),'-'); hold on;
    plot(xc(:,i),yc(:,i),'-'); hold on;
end
xp=zeros(20,1);
yp=zeros(20,1);
xp=0:0.05:1;
yp=0.8
yp2=0.3*cos(pi*xp)+0.5;
yp3=0.2;
yp4=-0.3*cos(pi*xp)+0.5;
plot(xp,yp,'or');
text(0.05,0.12,'F');
% FACE II
for i=1:2,
    plot(xc(i,:)+1,yc(i,:),'-'); hold on;
    plot(xc(:,i)+1,yc(:,i),'-'); hold on;
end
plot(xp+1,yp2,'or');
text(1.88,0.90,'E');
% FACE III
for i=1:2,
    plot(xc(i,:)+2,yc(i,:),'-'); hold on;
    plot(xc(:,i)+2,yc(:,i),'-'); hold on;
end
plot(xp+2,yp3,'or');
text(1+1.88,0.90,'B');
% FACE IV
for i=1:2,
    plot(xc(i,:)-1,yc(i,:),'-'); hold on;
    plot(xc(:,i)-1,yc(:,i),'-'); hold on;
end
plot(xp-1,yp4,'or');
text(0.88-1-0.8,0.90,'W');
% FACE V
for i=1:2,
    plot(xc(i,:),yc(i,:)+1,'-'); hold on;
    plot(xc(:,i),yc(:,i)+1,'-'); hold on;
end
% FACE VI
for i=1:2,
    plot(xc(i,:),yc(i,:)-1,'-'); hold on;
    plot(xc(:,i),yc(:,i)-1,'-'); hold on;
end
axis([-1.5 3.5 -1.5 2.5]);
title('NOTATION FOR THE 6 PATCHES OF THE CUBED-SPHERE');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
for i=1:nn,
  for j=1:nn,
    xi1(i,j)=xi(i);
    eta1(i,j)=eta(j);
  end
end
for j=1:nn,
    plot(xi1(:,j),eta1(:,j),'k'); hold on;
end
% PLOT VERTICAL ISOLINES
for i=1:nn,
    plot(xi1(i,:),eta1(i,:),'k'); hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure(6);
axis([-1 1 -1 1]);
amin=-pi/4;
amax=pi/4; 
na=16;
da=pi/(2*na);
t = amin:da:amax;
[xi2,eta2] = meshgrid(t);
set(gcf,'Color',[0.5,0.5,0.5])
subplot(2,2,1)
for j=1:na+1,
    plot(xi2(:,j),eta2(:,j),'k'); hold on;
end
size(xi2)
for i=1:na+1,
    plot(xi2(i,:),eta2(i,:),'k'); hold on;
end
nnn=13
%for i=1:na+1,
% '--ko'== black circle
% '--rs'== red square
plot(xi2(nnn,:),eta2(nnn,:),'--k>','LineWidth',1,...
                'MarkerEdgeColor',[0.2 0.2 0.2],...
                'MarkerFaceColor',[0.2 0.2 0.2],...
                'MarkerSize',6)
%plot(xi2(13,4),eta2(13,4));          
%end
%str1=strcat('p','/4');
%str=['p';'2p';'3p';'4p';'5p'];
%str=[str1;'2p';'3p';'4p';'5p'];
% set(gca,'xticklabel',str,'fontname','symbol')
%str=['-0.25p';'0';'0.25p'];
%set(gca,'xticklabel',str,'fontname','symbol')
%set(gca,'XTick',-pi/4:pi/4:pi/4)
%set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'})
set(gca,'XTick',-pi/4:pi/4:pi/4)
set(gca,'XTickLabel',{'-pi/4','0','pi/4'})
xlabel('-\pi/4 \leq \xi_F \leq \pi/4')
set(gca,'YTick',-pi/4:pi/4:pi/4)
set(gca,'YTickLabel',{'-\pi/4','0','\pi/4'})
ylabel('-\pi/4 \leq \eta_F \leq \pi/4')
title('FRONT PATCH');
text(-pi/2.8,pi/7.5,'\eta=\eta_0','FontSize',10)
axis equal
%%%%%%%%%%%%%%%%
subplot(2,2,2)
for j=1:na+1,
    plot(xi2(:,j),eta2(:,j),'k'); hold on;
end
size(xi2)
for i=1:na+1,
    plot(xi2(i,:),eta2(i,:),'k'); hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL CROSS EN XI/ETA SUR FACE EST
for i=1:na+1,
    beta3(i)=atan(-sin(xi(i))*tan(eta(nnn)));
    eta3(i)=atan(sqrt(1+tan(xi(i))^2)*tan(beta3(i)));
end
plot(xi,eta3,'--k>','LineWidth',1,...
                'MarkerEdgeColor',[0.2 0.2 0.2],...
                'MarkerFaceColor',[0.2 0.2 0.2],...
                'MarkerSize',6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5            

set(gca,'XTick',-pi/4:pi/4:pi/4)
set(gca,'XTickLabel',{'-pi/4','0','pi/4'})
xlabel('-\pi/4 \leq \xi_E \leq \pi/4')
set(gca,'YTick',-pi/4:pi/4:pi/4)
set(gca,'YTickLabel',{'-pi/4','0','pi/4'})
ylabel('-\pi/4 \leq \eta_E \leq \pi/4')
title('EAST PATCH');
axis equal
%%%%%%%%%%%%%%%%5
subplot(2,2,3)
for j=1:na+1,
    plot(xi2(:,j),eta2(:,j),'k'); hold on;
end
size(xi2)
for i=1:na+1,
    plot(xi2(i,:),eta2(i,:),'k'); hold on;
end
set(gca,'XTick',-pi/4:pi/4:pi/4)
set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'})
xlabel('-\pi/4 \leq \xi_B \leq \pi/4')
set(gca,'YTick',-pi/4:pi/4:pi/4)
set(gca,'YTickLabel',{'-\pi/4','0','\pi/4'})
ylabel('-\pi/4 \leq \eta_B \leq \pi/4')
axis equal
title('BACK PATCH');
plot(xi2(nnn,:),eta2(na+1-nnn,:),'--k>','LineWidth',1,...
                'MarkerEdgeColor',[0.2 0.2 0.2],...
                'MarkerFaceColor',[0.2 0.2 0.2],...
                'MarkerSize',6)
%%%%%%%%%%%%%%%%5
subplot(2,2,4)
for j=1:na+1,
    plot(xi2(:,j),eta2(:,j),'k'); hold on;
end
size(xi2)
for i=1:na+1,
    plot(xi2(i,:),eta2(i,:),'k'); hold on;
end
set(gca,'XTick',-pi/4:pi/4:pi/4)
set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'})
xlabel('-\pi/4 \leq \xi_W \leq \pi/4')
set(gca,'YTick',-pi/4:pi/4:pi/4)
set(gca,'YTickLabel',{'-\pi/4','0','\pi/4'})
ylabel('-\pi/4 \leq \eta_W \leq \pi/4')
axis equal
title('WEST PATCH');
%%%%%%%%%%%%%%%%
%set(findobj(gca,'Type','line','Color',[0 0 1]),...
 %   'Color','grey',...
  %  'LineWidth',0.5)
%annotation('arrow',x,y)
axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL CROSS EN XI/ETA SUR FACE EST
na+1-nnn
for i=1:na+1,
    beta4(i)=atan(-sin(xi(i))*tan(eta(na+1-nnn)));
    eta4(i)=atan(sqrt(1+tan(xi(i))^2)*tan(beta4(i)));
end
plot(xi,eta4,'--k>','LineWidth',1,...
                'MarkerEdgeColor',[0.2 0.2 0.2],...
                'MarkerFaceColor',[0.2 0.2 0.2],...
                'MarkerSize',6)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7- DESSIN DES 6 CARRES POUR ARTCLE CS1
% ET EXPOSE TALKS/CS.
% un carre
figure(7);
xc=zeros(2,2);yc=zeros(2,2);
for j=1:2,
        xc(:,j)=0:1;
        yc(j,:)=0:1;
end
% FACE I
for i=1:2,
    plot(xc(i,:),yc(i,:),'-'); hold on;
    plot(xc(:,i),yc(:,i),'-'); hold on;
end
xp=zeros(20,1);
yp=zeros(20,1);
xp=0:0.05:1;
yp=0.8*ones(21,1);
size(yp)
yp2=0.3*cos(pi*xp)+0.5;
yp3=0.2;
yp4=-0.3*cos(pi*xp)+0.5;
%plot(xp,yp,'or');
text(0.05,0.12,'Front');
size(xp)
size(yp)
plot(xp,yp,'kx'); hold on;
plot(xp,yp,'or');
% FACE II
for i=1:2,
    plot(xc(i,:)+1,yc(i,:),'-'); hold on;
    plot(xc(:,i)+1,yc(:,i),'-'); hold on;
end
plot(xp+1,yp2,'or');
text(1.88,0.90,'East');
xp5=zeros(20,1);
for i=1:21,
    xp5=(xp(i)+1)*ones(21,1);
    yp5=0:0.05:1;
    plot(xp5,yp5,'k-');hold on;
end
% FACE III
for i=1:2,
    plot(xc(i,:)+2,yc(i,:),'-'); hold on;
    plot(xc(:,i)+2,yc(:,i),'-'); hold on;
end
plot(xp+2,yp3,'or');
plot(xp+2,yp3,'kx');
text(1+1.88,0.90,'Back');
% FACE IV
for i=1:2,
    plot(xc(i,:)-1,yc(i,:),'-'); hold on;
    plot(xc(:,i)-1,yc(:,i),'-'); hold on;
end
plot(xp-1,yp4,'or');
text(0.88-1-0.8,0.90,'West');
for i=1:21,
    xp6=(xp(i)-1)*ones(21,1);
    yp6=0:0.05:1;
    plot(xp6,yp6,'k-');hold on;
end
% FACE V
for i=1:2,
    plot(xc(i,:),yc(i,:)+1,'-'); hold on;
    plot(xc(:,i),yc(:,i)+1,'-'); hold on;
end
% FACE VI
for i=1:2,
    plot(xc(i,:),yc(i,:)-1,'-'); hold on;
    plot(xc(:,i),yc(:,i)-1,'-'); hold on;
end
axis([-1.5 3.5 -1.5 2.5]);
title('SHAPE OF A GREAT CIRCLE OF THE NETWORK I_{alpha}');
% figure (8)
%% PLOT OF THE PROJECTION OF THE GRID ON FACE I ON THE CUBE.
%
figure(8);
yc=zeros(nn);
zc=zeros(nn);
vv=zeros(nn);
%
yc=tan(xi);
zc=tan(eta);
surf(yc,zc,vv);
colormap(gray);
title('PATCH I OF THE CUBED-SPHERE MAPPED ONTO THE CUBE');
set(gca,'XTick',-1:0.5:1);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1'});
xlabel('-1 \leq y \leq 1');
set(gca,'YTick',-1:0.5:1);
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'});
ylabel('-1 \leq z \leq 1');
view(2);
% text(-pi/2.8,pi/7.5,'\eta=\eta_0','FontSize',10)
axis equal;
% figure (9)
%% PLOT OF THE PROJECTION OF THE GRID ON FACE I ON THE CUBE.
%
figure(9);
title('PATCH FRONT, MAPPED GRID ON FACE FRONT OF THE CUBE'); 
yc=zeros(nn);
zc=zeros(nn);
vv=zeros(nn);
%
subplot(1,2,1);
surf(xi,eta,vv);
%colormap(gray);
colormap(jet);
set(gca,'XTick',-pi/4:pi/4:pi/4);
set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'});
xlabel('-\pi/4 \leq \xi_F \leq \pi/4');
set(gca,'YTick',-pi/4:pi/4:pi/4);
set(gca,'YTickLabel',{'-\pi/4','0','\pi/4'});
ylabel('-\pi/4 \leq \eta_F \leq \pi/4');
%title('PATCH FRONT, CARTESIAN GRID IN THE (\xi,\eta) SPACE'); 
view(2);
axis equal;
%%%
subplot(1,2,2);
yc=tan(xi);
zc=tan(eta);
surf(yc,zc,vv);
%colormap(gray);
colormap(jet);
set(gca,'XTick',-1:0.5:1);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1'});
xlabel('-1 \leq y \leq 1');
set(gca,'YTick',-1:0.5:1);
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'});
ylabel('-1 \leq z \leq 1');
%title('PATCH FRONT, MAPPED GRID ON FACE FRONT OF THE CUBE'); 
view(2);
% text(-pi/2.8,pi/7.5,'\eta=\eta_0','FontSize',10)
axis equal;
% 10- DESSIN DE LA TOPOLOGIE DE LA CUBED-SPHERE
figure(10);
xc=zeros(2,2);yc=zeros(2,2);
for j=1:2,
        xc(:,j)=0:1;
        yc(j,:)=0:1;
end
% FACE I
%for i=1:2,
plot(xc(1,:),yc(1,:),'k-');
%set(gca1,'Color','red');    
hold on;
%    set(gca,'Color','red');
plot(xc(:,1),yc(:,1),'k-'); 
hold on;
%%%
plot(xc(2,:),yc(2,:),'k-');
%set(gca1,'Color','red');    
hold on;
%    set(gca,'Color','red');
plot(xc(:,2),yc(:,2),'k-'); 
hold on;
%%
text(0.3,0.5,'FRONT');
%
% FACE II
%for i=1:2,
plot(xc(1,:)+1,yc(1,:),'-'); hold on;
%
plot(xc(:,1)+1,yc(:,1),'-','Color','red','Linewidth',4); hold on;
%
plot(xc(2,:)+1,yc(2,:)); hold on;
%
plot(xc(:,2)+1,yc(:,2),'-','Color','green','Linewidth',4); hold on;
%%%
%end
text(0.3+1,0.5,'EAST');
% FACE III
plot(xc(1,:)+2,yc(1,:),'-'); hold on;
plot(xc(:,1)+2,yc(:,1),'Color','blue','Linewidth',4); hold on;
plot(xc(2,:)+2,yc(2,:),'Color','black','Linewidth',4); hold on;
plot(xc(:,2)+2,yc(:,2),'Color','yellow','Linewidth',4); hold on;
text(0.3+2,0.5,'BACK');
% FACE IV
plot(xc(1,:)-1,yc(1,:),'Color','black','Linewidth',4); hold on;
plot(xc(:,1)-1,yc(:,1),'Color','magenta','Linewidth',4); hold on;
plot(xc(2,:)-1,yc(2,:),'-'); hold on;
plot(xc(:,2)-1,yc(:,2),'Color','cyan','Linewidth',4); hold on;
text(0.3-1,0.5,'WEST');
% FACE V
%for i=1:2,
plot(xc(1,:),yc(1,:)+1,'Color','cyan','Linewidth',4); hold on;
plot(xc(:,1),yc(:,1)+1,'-'); hold on;
plot(xc(2,:),yc(2,:)+1,'-','Color','green','Linewidth',4); hold on;
plot(xc(:,2),yc(:,2)+1,'-','Color','yellow','Linewidth',4); hold on;
%end
text(0.3,0.5+1,'NORTH');
% FACE VI
%for i=1:2,
plot(xc(1,:),yc(1,:)-1,'Color','magenta','Linewidth',4); hold on;
plot(xc(:,1),yc(:,1)-1,'-','Color','blue','Linewidth',4); hold on;
plot(xc(2,:),yc(2,:)-1,'-','Color','red','Linewidth',4); hold on;
plot(xc(:,2),yc(:,2)-1,'-'); hold on;
%end
text(0.3,0.5-1,'SOUTH');
axis([-1.5 3.5 -1.5 2.5]);
title('CUBED-SPHERE : TOPOLOGY AND PERIODICITY');
% ARROWS ARE ADDED BY HAND !!!!!!


            
