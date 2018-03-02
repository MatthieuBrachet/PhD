function [x,y] = equateur(ht_fI, ht_fII, ht_fIII, ht_fIV,ht_fV, ht_fVI)
% coupe le long de l'équateur
global nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI

pc=floor((nn+1)/2);
% panel IV
xx=x_fIV; yy=y_fIV; zz=z_fIV;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae1=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte1=ht_fIV(:,pc);

% panel I
xx=x_fI; yy=y_fI; zz=z_fI;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae2a=ll(1:pc);
hte2a=ht_fI(1:pc,pc);
lambdae2b=ll(pc+1:end);
hte2b=ht_fI(pc+1:end,pc);

% panel II
xx=x_fII; yy=y_fII; zz=z_fII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae3=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte3=ht_fII(:,pc);

% panel III
xx=x_fIII; yy=y_fIII; zz=z_fIII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae4=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte4=ht_fIII(:,pc);

% assemblage
x=[lambdae2b' lambdae3' lambdae4' lambdae1' lambdae2a'+2*pi];
y=[hte2b' hte3' hte4' hte1' hte2a'];
% figure(400)
% %plot(lambdae1,hte1,'r',lambdae2a,hte2a,'b',lambdae2b,hte2b,'c',lambdae3,hte3,'g',lambdae4,hte4,'y')
% plot(x,y)
% title('Equator')
% grid on
end

