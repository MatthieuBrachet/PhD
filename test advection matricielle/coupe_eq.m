function [x,f] = coupe_eq(funfI,funfII,funfIII,funfIV)
global nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
% [lambda, ~, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
% p=floor(size(lambda,1)/2);
% f=funfIV(:,p);
% x=lambda(:,p);

pc=floor((nn+1)/2);
% panel IV
xx=x_fIV; yy=y_fIV; zz=z_fIV;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae1=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte1=funfIV(:,pc);

% panel I
xx=x_fI; yy=y_fI; zz=z_fI;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae2a=ll(1:pc);
hte2a=funfI(1:pc,pc);
lambdae2b=ll(pc+1:end);
hte2b=funfI(pc+1:end,pc);

% panel II
xx=x_fII; yy=y_fII; zz=z_fII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae3=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte3=funfII(:,pc);

% panel III
xx=x_fIII; yy=y_fIII; zz=z_fIII;
[lambda, teta, ~]=cart2sph(xx,yy,zz);
ll=lambda(:,pc);
lambdae4=ll.*(ll>=0)+(ll+2*pi).*(ll<0);
hte4=funfIII(:,pc);

% assemblage
x=[lambdae3' lambdae4' lambdae1' lambdae2a'+2*pi];
f=[hte3' hte4' hte1' hte2a'];

x=[lambdae4(pc:end)' lambdae1' lambdae2a'+2*pi];
f=[hte4(pc:end)' hte1' hte2a'];

x=x(floor(end/2)+1:end);
f=f(floor(end/2)+1:end);
end