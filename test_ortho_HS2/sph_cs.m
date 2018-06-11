function [ hs ] = sph_cs( nhs, mhs )
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

%% x coord
% points interieurs (6)
xi_fI=reshape(x_fI(2:end-1,2:end-1),[],1);
xi_fII=reshape(x_fII(2:end-1,2:end-1),[],1);
xi_fIII=reshape(x_fIII(2:end-1,2:end-1),[],1);
xi_fIV=reshape(x_fIV(2:end-1,2:end-1),[],1);
xi_fV=reshape(x_fV(2:end-1,2:end-1),[],1);
xi_fVI=reshape(x_fVI(2:end-1,2:end-1),[],1);

xi=[xi_fI; xi_fII; xi_fIII; xi_fIV; xi_fV; xi_fVI];

% points de bords (12)
xb_1=reshape(x_fI(2:end-1,1),[],1);
xb_2=reshape(x_fII(2:end-1,1),[],1);
xb_3=reshape(x_fIII(2:end-1,1),[],1);
xb_4=reshape(x_fIV(2:end-1,1),[],1);

xb_5=reshape(x_fI(2:end-1,end),[],1);
xb_6=reshape(x_fII(2:end-1,end),[],1);
xb_7=reshape(x_fIII(2:end-1,end),[],1);
xb_8=reshape(x_fIV(2:end-1,end),[],1);

xb_9=reshape(x_fI(1,2:end-1),[],1);
xb_10=reshape(x_fII(1,2:end-1),[],1);
xb_11=reshape(x_fIII(1,2:end-1),[],1);
xb_12=reshape(x_fIV(1,2:end-1),[],1);

xb=[xb_1;xb_2;xb_3;xb_4;xb_5;xb_6;xb_7;xb_8;xb_9;xb_10;xb_11;xb_12];

% coins (8)
xc_1=x_fI(1,1);
xc_2=x_fI(1,end);
xc_3=x_fI(end,1);
xc_4=x_fI(end,end);

xc_5=x_fIII(1,1);
xc_6=x_fIII(1,end);
xc_7=x_fIII(end,1);
xc_8=x_fIII(end,end);

xc=[xc_1;xc_2;xc_3;xc_4;xc_5;xc_6;xc_7;xc_8];

% points
XX=[xi; xb; xc];

%% y coord
% points interieurs (6)
yi_fI=reshape(y_fI(2:end-1,2:end-1),[],1);
yi_fII=reshape(y_fII(2:end-1,2:end-1),[],1);
yi_fIII=reshape(y_fIII(2:end-1,2:end-1),[],1);
yi_fIV=reshape(y_fIV(2:end-1,2:end-1),[],1);
yi_fV=reshape(y_fV(2:end-1,2:end-1),[],1);
yi_fVI=reshape(y_fVI(2:end-1,2:end-1),[],1);

yi=[yi_fI; yi_fII; yi_fIII; yi_fIV; yi_fV; yi_fVI];

% points de bords (12)
yb_1=reshape(y_fI(2:end-1,1),[],1);
yb_2=reshape(y_fII(2:end-1,1),[],1);
yb_3=reshape(y_fIII(2:end-1,1),[],1);
yb_4=reshape(y_fIV(2:end-1,1),[],1);

yb_5=reshape(y_fI(2:end-1,end),[],1);
yb_6=reshape(y_fII(2:end-1,end),[],1);
yb_7=reshape(y_fIII(2:end-1,end),[],1);
yb_8=reshape(y_fIV(2:end-1,end),[],1);

yb_9=reshape(y_fI(1,2:end-1),[],1);
yb_10=reshape(y_fII(1,2:end-1),[],1);
yb_11=reshape(y_fIII(1,2:end-1),[],1);
yb_12=reshape(y_fIV(1,2:end-1),[],1);

yb=[yb_1;yb_2;yb_3;yb_4;yb_5;yb_6;yb_7;yb_8;yb_9;yb_10;yb_11;yb_12];

% coins (8)
yc_1=y_fI(1,1);
yc_2=y_fI(1,end);
yc_3=y_fI(end,1);
yc_4=y_fI(end,end);

yc_5=y_fIII(1,1);
yc_6=y_fIII(1,end);
yc_7=y_fIII(end,1);
yc_8=y_fIII(end,end);

yc=[yc_1;yc_2;yc_3;yc_4;yc_5;yc_6;yc_7;yc_8];

% points
YY=[yi; yb; yc];

%% z coord
% points interieurs (6)
zi_fI=reshape(z_fI(2:end-1,2:end-1),[],1);
zi_fII=reshape(z_fII(2:end-1,2:end-1),[],1);
zi_fIII=reshape(z_fIII(2:end-1,2:end-1),[],1);
zi_fIV=reshape(z_fIV(2:end-1,2:end-1),[],1);
zi_fV=reshape(z_fV(2:end-1,2:end-1),[],1);
zi_fVI=reshape(z_fVI(2:end-1,2:end-1),[],1);

zi=[zi_fI; zi_fII; zi_fIII; zi_fIV; zi_fV; zi_fVI];

% points de bords (12)
zb_1=reshape(z_fI(2:end-1,1),[],1);
zb_2=reshape(z_fII(2:end-1,1),[],1);
zb_3=reshape(z_fIII(2:end-1,1),[],1);
zb_4=reshape(z_fIV(2:end-1,1),[],1);

zb_5=reshape(z_fI(2:end-1,end),[],1);
zb_6=reshape(z_fII(2:end-1,end),[],1);
zb_7=reshape(z_fIII(2:end-1,end),[],1);
zb_8=reshape(z_fIV(2:end-1,end),[],1);

zb_9=reshape(z_fI(1,2:end-1),[],1);
zb_10=reshape(z_fII(1,2:end-1),[],1);
zb_11=reshape(z_fIII(1,2:end-1),[],1);
zb_12=reshape(z_fIV(1,2:end-1),[],1);

zb=[zb_1;zb_2;zb_3;zb_4;zb_5;zb_6;zb_7;zb_8;zb_9;zb_10;zb_11;zb_12];

% coins (8)
zc_1=z_fI(1,1);
zc_2=z_fI(1,end);
zc_3=z_fI(end,1);
zc_4=z_fI(end,end);

zc_5=z_fIII(1,1);
zc_6=z_fIII(1,end);
zc_7=z_fIII(end,1);
zc_8=z_fIII(end,end);

zc=[zc_1;zc_2;zc_3;zc_4;zc_5;zc_6;zc_7;zc_8];

% points
ZZ=[zi; zb; zc];

%% Calcul de l'harmonique
hs = sph( nhs,mhs,XX,YY,ZZ );
% 
% figure(100)
% X=reshape(XX,[],2);
% Y=reshape(YY,[],2);
% Z=reshape(ZZ,[],2);
% h=reshape(hs,[],2);
% surf(X,Y,Z,real(h));
% shading interp

end

