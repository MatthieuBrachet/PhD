%%% INTERPOLATION WITH SPHs ON THE WHOLE CUBED-SPHERE
%%% => NO PROPERTIES OF SYMMETRY FOR W_{i,j}
%%% => NOT THE SAME FOR THE 6 FACES

clear all;
global n nn;
global radius;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=32;
nhs_max=ceil(sqrt(6*N^2+2)-1);

make_cs_grid(N);

s=N+1;
xx_fI=x_fI; xx_fII=x_fII'; xx_fIII=x_fIII; xx_fIV=x_fIV'; xx_fV=-x_fV'; xx_fVI=-x_fVI';
yy_fI=y_fI; yy_fII=y_fII'; yy_fIII=y_fIII; yy_fIV=y_fIV'; yy_fV=-y_fV'; yy_fVI=-y_fVI';
zz_fI=z_fI; zz_fII=z_fII'; zz_fIII=z_fIII; zz_fIV=z_fIV'; zz_fV=-z_fV'; zz_fVI=-z_fVI';

xx=[reshape(xx_fI(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(xx_fII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(xx_fIII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(xx_fIV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(xx_fV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(xx_fVI(2:s-1,2:s-1),[(s-2)^2,1])];
x1=xx_fI(:,1); x2=xx_fI(:,s); x3=xx_fIII(:,1); x4=xx_fIII(:,s);
x5=xx_fV(s,:); x6=xx_fV(:,s); x7=xx_fV(1,:); x8=xx_fV(:,1);
x9=xx_fVI(1,:); x10=xx_fVI(:,s); x11=xx_fVI(s,:); x12=xx_fVI(:,1);
x_aretes1=[x1;x2;x3;x4];
x_aretes2=[x5(2:s-1)';x6(2:s-1);x7(2:s-1)';x8(2:s-1); ...
           x9(2:s-1)';x10(2:s-1);x11(2:s-1)';x12(2:s-1)];
xx=[xx;x_aretes1;x_aretes2];

yy=[reshape(yy_fI(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(yy_fII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(yy_fIII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(yy_fIV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(yy_fV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(yy_fVI(2:s-1,2:s-1),[(s-2)^2,1])];
y1=yy_fI(:,1); y2=yy_fI(:,s); y3=yy_fIII(:,1); y4=yy_fIII(:,s);
y5=yy_fV(s,:); y6=yy_fV(:,s); y7=yy_fV(1,:); y8=yy_fV(:,1);
y9=yy_fVI(1,:); y10=yy_fVI(:,s); y11=yy_fVI(s,:); y12=yy_fVI(:,1);
y_aretes1=[y1;y2;y3;y4];
y_aretes2=[y5(2:s-1)';y6(2:s-1);y7(2:s-1)';y8(2:s-1); ...
           y9(2:s-1)';y10(2:s-1);y11(2:s-1)';y12(2:s-1)];
yy=[yy;y_aretes1;y_aretes2];

zz=[reshape(zz_fI(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(zz_fII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(zz_fIII(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(zz_fIV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(zz_fV(2:s-1,2:s-1),[(s-2)^2,1]); ...
    reshape(zz_fVI(2:s-1,2:s-1),[(s-2)^2,1])];
z1=zz_fI(:,1); z2=zz_fI(:,s); z3=zz_fIII(:,1); z4=zz_fIII(:,s);
z5=zz_fV(s,:); z6=zz_fV(:,s); z7=zz_fV(1,:); z8=zz_fV(:,1);
z9=zz_fVI(1,:); z10=zz_fVI(:,s); z11=zz_fVI(s,:); z12=zz_fVI(:,1);
z_aretes1=[z1;z2;z3;z4];
z_aretes2=[z5(2:s-1)';z6(2:s-1);z7(2:s-1)';z8(2:s-1); ...
           z9(2:s-1)';z10(2:s-1);z11(2:s-1)';z12(2:s-1)];
zz=[zz;z_aretes1;z_aretes2];

A=[];
for nhs=0:1:nhs_max,
        for mhs=0:1:nhs, 
            nhs
            fun=sph(nhs,mhs,xx,yy,zz);
            A=[A;fun'];
        end
end  
%k=min(size(A));
k=ceil(size(A,1)/2);
b=zeros(k,1);
b(1)=sqrt(4*pi);
w=pinv(A(1:k,:))*b;

disp 'erreur fonction = 1'
abs(4*pi-sum(w))
disp 'erreur f1'
fun=fun_f1(xx,yy,zz);
abs(216*pi/35-fun'*w)
disp 'erreur f2'
fun=fun_f2(xx,yy,zz);
abs(6.6961822200736179523-fun'*w)
disp 'erreur f3'
fun=fun_f3(xx,yy,zz);
abs(4*pi/9-fun'*w)
disp 'erreur f4'
fun=fun_f4(xx,yy,zz);
abs(4*pi/9-fun'*w)

 
