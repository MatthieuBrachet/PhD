clear all;
%ERROR PICTURES FOR THE CASE kp=5 NUMBER OF BRANCHES
% LOWER GRAPHICS: LEAST SQUARE OF THE ERROR FOR PSI
% UPPER GRAPHIS: LEAST SQUARE OF THE ERROR FOR PSI_X
% hh=zeros(4,1);
% nn=[9;17;33;65;129];
fid=fopen('flower5b.txt','r');
frewind(fid);
%%
nline=6;
for i=1:nline,
    st=fgetl(fid);
    line1=sscanf(st,'%f',3);
    nn(i)=line1(1);
    ee(i)=line1(2);
    eex(i)=line1(3);
end
%
% nn=[9;17;33;65];
a=0; b=1;
hh=(b-a)./(nn-1);
lhh=log10(hh);
% ee=[1.7561e-5;3.9014e-7;1.7420e-8;9.9940e-10];
% eex=[4.4761e-5;3.5979e-6;9.5902e-7;9.1149e-8];
lee=log10(ee);
leex=log10(eex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11); % PLOT OF CONVERGENCE FOR PSI AND PSI_X IN MAX NORMS
% hdl1=plot(lhh,lee,'or',lhh,lee,'-b');
hl1=plot(lhh,lee,'or');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');
set(hl1,'MarkerFaceColor','red');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ht=text('Position',[-1.2,-6.8,0],'String','slope for \psi: 4.68');
% set(ht,'FontSize',12);
% set(ht,'FontWeight','demi');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ht=text('Position',[-1.8,-5.3,0],'String','slope for \psi_x: 2.87');
% set(ht,'FontSize',12);
% set(ht,'FontWeight','demi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
[a1,b1]=polyfit(lhh,lee,1);
tt=linspace(min(lhh),max(lhh),1000);
lsq1=a1(2)+a1(1)*tt;
hlf1=plot(tt,lsq1,'-b');
set(hlf1,'LineWidth',2.0);
hold on;
%%%%%%
hdl2=plot(lhh,leex,'or');
set(hdl2,'LineWidth',2.0);
set(hdl2,'MarkerSize',10);
set(hdl2,'MarkerEdgeColor','k');
set(hdl2,'MarkerFaceColor','green');
hold on;
[a2,b2]=polyfit(lhh,leex,1);
%tt=linspace(min(lhh),max(lhh),1000);
lsq2=a2(2)+a2(1)*tt;
hlf2=plot(tt,lsq2,'-b');
set(hlf2,'LineWidth',2.0);
hold on;
axis([-2 -0.5 -10 -3]);grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha=gca;
xa=get(gca,'Xlabel'); %%% POURQUOI "GCA" ET PAS HA ??
set(ha,'Xgrid','on');
set(xa,'String','Log10(h)');
set (xa,'FontName','Calibri');
set(xa,'FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ha=gca;
ya=get(gca,'Ylabel'); %%% POURQUOI "GCA" ET PAS HA ??
set(ha,'Ygrid','on');
set(ya,'String','Log10(error_{max})');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TO AJUST AT THE END %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
a1
%%% WRITE VALUE OF A1(1) IN THE NEXT LINE= ESTIMATED CONVERGENCE RATE FOR
%%% PSI. ALSO ADJUST POSITION OF THE TEXT
ht=text('Position',[-1.1,-6.8,0],'String','slope for \psi: 3.72');
set(ht,'FontSize',12);
set(ht,'FontWeight','demi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a2
%%% WRITE VALUE OF A2(1) IN THE NEXT LINE= ESTIMATED CONVERGENCE RATE FOR
%%% PSI_X. ALSO ADJUST POSITION OF THE TEXT
ht=text('Position',[-1.7,-5.3,0],'String','slope for \psi_x: 3.48');
set(ht,'FontSize',12);
set(ht,'FontWeight','demi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha=gca;
xt=get(gca,'Title');
title('Convergence rate: star with 5 branches');
set (xt,'FontSize',12);


