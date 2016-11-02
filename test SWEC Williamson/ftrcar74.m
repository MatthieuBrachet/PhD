function [ hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI, vf_fI, vf_fII, vf_fIII, vf_fIV, vf_fV, vf_fVI ]...
    = ftrcar74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
%% *** caracteristic filter
% caracteristic filter for the SWE problem : 
%
%   dv/dt + (v*nabla)*v + f * k vect v + gp*grad(h) = 0
%   dh/dt + div(h*v) = 0
%
% created the 26/10/2016 by Matthieu Brachet.

global gp
global nn n
% -----------------------------------
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gdxi_I gdxi_II gdxi_III gdxi_IV gdxi_V gdxi_VI;
global gdeta_I gdeta_II gdeta_III gdeta_IV gdeta_V gdeta_VI;
% -----------------------------------
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI

%% *** caracteristics variables *******************************************
[hs_fI] = relief(x_fI,y_fI,z_fI);
[hs_fII] = relief(x_fII,y_fII,z_fII);
[hs_fIII] = relief(x_fIII,y_fIII,z_fIII);
[hs_fIV] = relief(x_fIV,y_fIV,z_fIV);
[hs_fV] = relief(x_fV,y_fV,z_fV);
[hs_fVI] = relief(x_fVI,y_fVI,z_fVI);


% Face I
for i=1:nn
    for j=1:nn
        h=ht_fI(i,j)-hs_fI(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fI(i,j,1:3),gxi_I(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fI(i,j,1:3),geta_I(i,j,1:3));
        X1(3)=ht_fI(i,j);
        Y_fI(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end

% Face II
for i=1:nn
    for j=1:nn
        h=ht_fII(i,j)-hs_fII(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fII(i,j,1:3),gxi_I(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fII(i,j,1:3),geta_I(i,j,1:3));
        X1(3)=ht_fII(i,j);
        Y_fII(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end

% Face III
for i=1:nn
    for j=1:nn
        h=ht_fIII(i,j)-hs_fIII(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fIII(i,j,1:3),gxi_III(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fIII(i,j,1:3),geta_III(i,j,1:3));
        X1(3)=ht_fIII(i,j);
        Y_fIII(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end

% Face IV
for i=1:nn
    for j=1:nn
        h=ht_fIV(i,j)-hs_fIV(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fIV(i,j,1:3),gxi_IV(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fIV(i,j,1:3),geta_IV(i,j,1:3));
        X1(3)=ht_fIV(i,j);
        Y_fIV(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end

% Face V
for i=1:nn
    for j=1:nn
        h=ht_fV(i,j)-hs_fV(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fV(i,j,1:3),gxi_V(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fV(i,j,1:3),geta_V(i,j,1:3));
        X1(3)=ht_fV(i,j);
        Y_fV(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end

% Face VI
for i=1:nn
    for j=1:nn
        h=ht_fVI(i,j)-hs_fVI(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        X1(1)=dot(vt_fVI(i,j,1:3),gxi_VI(i,j,1:3));  %%% *** changement a vérifier !!!
        X1(2)=dot(vt_fVI(i,j,1:3),geta_VI(i,j,1:3));
        X1(3)=ht_fVI(i,j);
        Y_fVI(i,j,1:3)=V_xi\(V_eta\transpose(X1));
    end
end



%% *** Filter *************************************************************
[Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1)]=ftr74(Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1),n,nn);
[Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2)]=ftr74(Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2),n,nn);
[Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3)]=ftr74(Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3),n,nn);

%% *** back to classical variables ****************************************
% Face I
for i=1:nn
    for j=1:nn
        h=ht_fI(i,j)-hs_fI(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fI(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fI(i,j,1:3)=X(1).*gxi_I(i,j,1:3)+X(2)*geta_I(i,j,1:3);
        hf_fI(i,j)=X(3);
    end
end

% Face II
for i=1:nn
    for j=1:nn
        h=ht_fII(i,j)-hs_fII(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fII(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fII(i,j,1:3)=X(1).*gxi_II(i,j,1:3)+X(2)*geta_II(i,j,1:3);
        hf_fII(i,j)=X(3);
    end
end

% Face III
for i=1:nn
    for j=1:nn
        h=ht_fIII(i,j)-hs_fIII(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fIII(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fIII(i,j,1:3)=X(1).*gxi_III(i,j,1:3)+X(2)*geta_III(i,j,1:3);
        hf_fIII(i,j)=X(3);
    end
end

% Face IV
for i=1:nn
    for j=1:nn
        h=ht_fIV(i,j)-hs_fIV(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fIV(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fIV(i,j,1:3)=X(1).*gxi_IV(i,j,1:3)+X(2)*geta_IV(i,j,1:3);
        hf_fIV(i,j)=X(3);
    end
end

% Face V
for i=1:nn
    for j=1:nn
        h=ht_fV(i,j)-hs_fV(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fV(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fV(i,j,1:3)=X(1).*gxi_V(i,j,1:3)+X(2)*geta_V(i,j,1:3);
        hf_fV(i,j)=X(3);
    end
end

% Face VI
for i=1:nn
    for j=1:nn
        h=ht_fVI(i,j)-hs_fVI(i,j);
        
        V11=0;
        V12=sqrt(gp*h);
        V13=-sqrt(gp*h);
        V21=h;
        V22=0;
        V23=0;
        V31=0;
        V32=h;
        V33=h;
        V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        V11=h;
        V12=0;
        V13=0;
        V21=0;
        V22=sqrt(h*gp);
        V23=-sqrt(h*gp);
        V31=0;
        V32=h;
        V33=h;
        V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fVI(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fVI(i,j,1:3)=X(1).*gxi_VI(i,j,1:3)+X(2)*geta_VI(i,j,1:3);
        hf_fVI(i,j)=X(3);
    end
end




end
