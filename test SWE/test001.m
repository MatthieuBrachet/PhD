mfunfI=vt_fI;
mfunfII=vt_fII;
mfunfIII=vt_fIII;
mfunfIV=vt_fIV;
mfunfV=vt_fV;
mfunfVI=vt_fVI;

global na nb;
global radius;
global alfa beta;
global betacr;
global alfa1;
global p_div k_div;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;

%% ************************************************************************ FRONT and BOTTOM on XI
xxtIa_I=zeros(nn,nn);xxtIa_II=zeros(nn,nn);
xxtIa_III=zeros(nn,nn);xxtIa_IV=zeros(nn,nn);
yytIa_I=zeros(nn,nn);yytIa_II=zeros(nn,nn);
yytIa_III=zeros(nn,nn);yytIa_IV=zeros(nn,nn);

% FACE I
for i=1:nn,
    for j=1:nn,
        xxtIa_I(i,j)=y_fI(i,j)/x_fI(i,j);
        yytIa_I(i,j)=z_fI(i,j)/x_fI(i,j);
    end
end
% FACE II
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fII(i,j));
        xwk2=fun7(x_fII(i,j));
        xxtIa_II(i,j)=y_fII(i,j)*xwk1;
        yytIa_II(i,j)=z_fII(i,j)*xwk2;
    end
end
% FACE III
for i=1:nn,
    for j=1:nn,
        xxtIa_III(i,j)=y_fIII(i,j)/x_fIII(i,j);
        yytIa_III(i,j)=-z_fIII(i,j)/x_fIII(i,j);
    end
end
% FACE IV
for i=1:nn,
    for j=1:nn,
        xwk1=fun5(x_fII(i,j)); 
        xwk2=fun7(x_fII(i,j)); 
        xxtIa_IV(i,j)=y_fIV(i,j)*xwk1;
        yytIa_IV(i,j)=z_fIV(i,j)*xwk2;
    end
end

deltatIa_I=(1+xxtIa_I.^2+yytIa_I.^2); 
deltatIa_II=(1+xxtIa_II.^2+yytIa_II.^2); 
deltatIa_III=(1+xxtIa_III.^2+yytIa_III.^2); 
deltatIa_IV=(1+xxtIa_IV.^2+yytIa_IV.^2);
deltabtIa_I=sqrt(deltatIa_I);
deltabtIa_II=sqrt(deltatIa_II);
deltabtIa_III=sqrt(deltatIa_III);
deltabtIa_IV=sqrt(deltatIa_IV);
gtIa_I=(radius^2)*(1+xxtIa_I.^2).*(1+yytIa_I.^2)./(deltabtIa_I.^3);
gtIa_II=(radius^2)*(1+xxtIa_II.^2).*(1+yytIa_II.^2)./(deltabtIa_II.^3);
gtIa_III=(radius^2)*(1+xxtIa_III.^2).*(1+yytIa_III.^2)./(deltabtIa_III.^3);
gtIa_IV=(radius^2)*(1+xxtIa_IV.^2).*(1+yytIa_IV.^2)./(deltabtIa_IV.^3);

gxiIa_II=zeros(nn,nn,3); gxiIa_IV=zeros(nn,nn,3);

% - FACE I -
gxiIa_I=gxi_I;
% - FACE II -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fII(i,j));
      xwk2=fun7(x_fII(i,j));
      gxiIa_II(i,j,1)= -y_fII(i,j)*xwk1*xwk1;
      gxiIa_II(i,j,2)= xwk1;
      gxiIa_II(i,j,3)= 0;
      
      gxiIa_II(i,j,1:3)=gxiIa_II(i,j,1:3)/(1+(y_fII(i,j)*xwk2)^2);
    end
end
% - FACE III -
gxiIa_III=gxi_III;
% - FACE IV -
for i=1:nn,
    for j=1:nn,
      xwk1=fun5(x_fIV(i,j));
      xwk2=fun7(x_fIV(i,j));
      gxiIa_IV(i,j,1)= -y_fIV(i,j)*xwk1*xwk1;
      gxiIa_IV(i,j,2)= xwk1;
      gxiIa_IV(i,j,3)= 0;
      %
      gxiIa_IV(i,j,1:3)=gxiIa_IV(i,j,1:3)/(1+(y_fIV(i,j)*xwk2)^2);
    end
end

funfI=zeros(nn,nn);funfII=zeros(nn,nn);funfIII=zeros(nn,nn);
funfIV=zeros(nn,nn);funfV=zeros(nn,nn);funfVI=zeros(nn,nn);
%
for i=1:nn,
    for j=1:nn,
      v1=mfunfI(i,j,:);
      v2=gxiIa_I(i,j,:);
      funfI(i,j,1:3)=cross(v1,v2);
      funfI(i,j,1:3)=funfI(i,j,1:3)*gtIa_I(i,j);
     end
end

