function [ lambdae, hte ] = equateur(ht_fI, ht_fII, ht_fIII, ht_fIV,ht_fV, ht_fVI)
% coupe le long de l'équateur
global nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI

pc=floor((nn+1)/2);
% panel IV
xx=x_fIV; yy=y_fIV; zz=z_fIV;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
lambdae1=lambda(:,pc);
hte1=ht_fIV(:,pc);

% panel I
xx=x_fI; yy=y_fI; zz=z_fI;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
lambdae2=lambda(:,pc);
hte2=ht_fI(:,pc);

% panel II
xx=x_fII; yy=y_fII; zz=z_fII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
lambdae3=lambda(:,pc);
hte3=ht_fII(:,pc);

% panel III
xx=x_fIII; yy=y_fIII; zz=z_fIII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
lambdae4=lambda(:,pc);
hte4=ht_fIII(:,pc);

% assemblage
lambdae=[lambdae1; lambdae2; lambdae3];
hte=[hte1; hte2; hte3];
end

