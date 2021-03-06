function[eq_fI, eq_fII, eq_fIII, eq_fIV, eq_fV, eq_fVI]=eq_cons101(ht_fI, ht_fII,...
    ht_fIII, ht_fIV, ht_fV, ht_fVI)
% conservation equation
global test
global n nn
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI
global radius

if test == 0
    %% test 3 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    %% Panel I
    xx=x_fI; yy=y_fI; zz=z_fI;
    u=ht_fI;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_I(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_I(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_I(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel II
    xx=x_fII; yy=y_fII; zz=z_fII;
    u=ht_fII;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_II(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_II(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_II(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel III
    xx=x_fIII; yy=y_fIII; zz=z_fIII;
    u=ht_fIII;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_III(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_III(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_III(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel IV
    xx=x_fIV; yy=y_fIV; zz=z_fIV;
    u=ht_fIV;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_IV(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_IV(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_IV(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel V
    xx=x_fV; yy=y_fV; zz=z_fV;
    u=ht_fV;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_V(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_V(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_V(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel VI
    xx=x_fVI; yy=y_fVI; zz=z_fVI;
    u=ht_fVI;
    f1=.5*u.^2;
    f2=.5*u.^2;
    f3=.5*u.^2;
    F_VI(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_VI(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_VI(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Assemblage
    [div_I,div_II,div_III,div_IV,div_V,div_VI]=...
            div101(F_I,F_II,F_III,F_IV,F_V,F_VI,n,nn);
        
    eq_fI=-div_I;
    eq_fII=-div_II;
    eq_fIII=-div_III;
    eq_fIV=-div_IV;
    eq_fV=-div_V;
    eq_fVI=-div_VI;
    
elseif test == 1
    %% test 1 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    %% Panel I
    xx=x_fI; yy=y_fI; zz=z_fI;
    u=ht_fI;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_I(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_I(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_I(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel II
    xx=x_fII; yy=y_fII; zz=z_fII;
    u=ht_fII;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_II(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_II(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_II(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel III
    xx=x_fIII; yy=y_fIII; zz=z_fIII;
    u=ht_fIII;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_III(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_III(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_III(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel IV
    xx=x_fIV; yy=y_fIV; zz=z_fIV;
    u=ht_fIV;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_IV(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_IV(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_IV(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel V
    xx=x_fV; yy=y_fV; zz=z_fV;
    u=ht_fV;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_V(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_V(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_V(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel VI
    xx=x_fVI; yy=y_fVI; zz=z_fVI;
    u=ht_fVI;
    f1=zeros(size(u));
    f2=zeros(size(u));
    f3=-pi*u.^2;
    F_VI(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_VI(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_VI(:,:,3)=(xx.*f2-yy.*f1)./radius;

    %% Assemblage
    [div_I,div_II,div_III,div_IV,div_V,div_VI]=...
            div101(F_I,F_II,F_III,F_IV,F_V,F_VI,n,nn);
        
    eq_fI=-div_I;
    eq_fII=-div_II;
    eq_fIII=-div_III;
    eq_fIV=-div_IV;
    eq_fV=-div_V;
    eq_fVI=-div_VI;
    
    
    
    
    elseif test == 2
    %% test 3 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    %% Panel I
    x=x_fI; y=y_fI; z=z_fI;
    u=ht_fI;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_I(:,:,1)=zeros(size(x));
    F_I(:,:,2)=z.*dh;
    F_I(:,:,3)=-y.*dh;
    
    %% Panel II
    x=x_fII; y=y_fII; z=z_fII;
    u=ht_fII;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_II(:,:,1)=zeros(size(x));
    F_II(:,:,2)=z.*dh;
    F_II(:,:,3)=-y.*dh;

    %% Panel III
    x=x_fIII; y=y_fIII; z=z_fIII;
    u=ht_fIII;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_III(:,:,1)=zeros(size(x));
    F_III(:,:,2)=z.*dh;
    F_III(:,:,3)=-y.*dh;
    
    %% Panel IV
    x=x_fIV; y=y_fIV; z=z_fIV;
    u=ht_fIV;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_IV(:,:,1)=zeros(size(x));
    F_IV(:,:,2)=z.*dh;
    F_IV(:,:,3)=-y.*dh;
    
    %% Panel V
    x=x_fV; y=y_fV; z=z_fV;
    u=ht_fV;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_V(:,:,1)=zeros(size(x));
    F_V(:,:,2)=z.*dh;
    F_V(:,:,3)=-y.*dh;
    
    %% Panel VI
    x=x_fVI; y=y_fVI; z=z_fVI;
    u=ht_fVI;
    dh=(fun1(x)+x.*dfun1(x)).*.5.*u.^2;
    F_VI(:,:,1)=zeros(size(x));
    F_VI(:,:,2)=z.*dh;
    F_VI(:,:,3)=-y.*dh;
    
    %% Assemblage
    [div_I,div_II,div_III,div_IV,div_V,div_VI]=...
            div101(F_I,F_II,F_III,F_IV,F_V,F_VI,n,nn);
        
    eq_fI=-div_I;
    eq_fII=-div_II;
    eq_fIII=-div_III;
    eq_fIV=-div_IV;
    eq_fV=-div_V;
    eq_fVI=-div_VI;
    
elseif test == 3
    %% test 3 of M. Ben-Artzi, J. Falcovitz and P. G. Lefloch
    %% Panel I
    xx=x_fI; yy=y_fI; zz=z_fI;
    f1=.5.*ht_fI.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_I(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_I(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_I(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel II
    xx=x_fII; yy=y_fII; zz=z_fII;
    f1=.5.*ht_fII.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_II(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_II(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_II(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel III
    xx=x_fIII; yy=y_fIII; zz=z_fIII;
    f1=.5.*ht_fIII.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_III(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_III(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_III(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel IV
    xx=x_fIV; yy=y_fIV; zz=z_fIV;
    f1=.5.*ht_fIV.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_IV(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_IV(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_IV(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel V
    xx=x_fV; yy=y_fV; zz=z_fV;
    f1=.5.*ht_fV.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_V(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_V(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_V(:,:,3)=(xx.*f2-yy.*f1)./radius;
    
    %% Panel VI
    xx=x_fVI; yy=y_fVI; zz=z_fVI;
    f1=.5.*ht_fVI.^2;
    f2=zeros(size(xx));
    f3=zeros(size(xx));
    F_VI(:,:,1)=(yy.*f3-zz.*f2)./radius;
    F_VI(:,:,2)=(zz.*f1-xx.*f3)./radius;
    F_VI(:,:,3)=(xx.*f2-yy.*f1)./radius;

    %% Assemblage
    [div_I,div_II,div_III,div_IV,div_V,div_VI]=...
            div101(F_I,F_II,F_III,F_IV,F_V,F_VI,n,nn);
        
    eq_fI=-div_I;
    eq_fII=-div_II;
    eq_fIII=-div_III;
    eq_fIV=-div_IV;
    eq_fV=-div_V;
    eq_fVI=-div_VI;
end
end

