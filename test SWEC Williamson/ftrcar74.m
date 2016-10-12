function [ hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI, vf_fI, vf_fII, vf_fIII, vf_fIV, vf_fV, vf_fVI ]...
    = ftrcar74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
%% *** caracteristic filter
% caracteristic filter for the L=SWE problem : 
%
%   dv/dt + div(0.5*norm(v) + g h ) + (f+vort(v))*k vect v = 0
%   dhs/dt + div(hs*v) = 0
%
% with hs=h-h_surf, h_surf is the relief function.
% created the 12/10/2016 by Matthieu Brachet.

global nn n
global G11_fI G12_fI G22_fI
global G11_fII G12_fII G22_fII
global G11_fIII G12_fIII G22_fIII
global G11_fIV G12_fIV G22_fIV
global G11_fV G12_fV G22_fV
global G11_fVI G12_fVI G22_fVI
% ----------------------------------
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gdxi_I gdxi_II gdxi_III gdxi_IV gdxi_V gdxi_VI;
global gdeta_I gdeta_II gdeta_III gdeta_IV gdeta_V gdeta_VI;
% -----------------------------------
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
% -----------------------------------
global gp

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
        uxi=dot(vt_fI(i,j,1:3),gdxi_I(i,j,1:3));
        ueta=dot(vt_fI(i,j,1:3),gdeta_I(i,j,1:3));
        G11=G11_fI(i,j);
        G12=G12_fI(i,j);
        G22=G22_fI(i,j);
        h=ht_fI(i,j);
        hs=hs_fI(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fI(i,j,1:3),gdxi_I(i,j,1:3));
        X(2)=dot(vt_fI(i,j,1:3),gdeta_I(i,j,1:3));
        X(3)=ht_fI(i,j);
        Y_fI(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face II
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fII(i,j,1:3),gdxi_II(i,j,1:3));
        ueta=dot(vt_fII(i,j,1:3),gdeta_II(i,j,1:3));
        G11=G11_fII(i,j);
        G12=G12_fII(i,j);
        G22=G22_fII(i,j);
        h=ht_fII(i,j);
        hs=hs_fII(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fII(i,j,1:3),gdxi_II(i,j,1:3));
        X(2)=dot(vt_fII(i,j,1:3),gdeta_II(i,j,1:3));
        X(3)=ht_fII(i,j);
        Y_fII(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face III
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fIII(i,j,1:3),gdxi_III(i,j,1:3));
        ueta=dot(vt_fIII(i,j,1:3),gdeta_III(i,j,1:3));
        G11=G11_fIII(i,j);
        G12=G12_fIII(i,j);
        G22=G22_fIII(i,j);
        h=ht_fIII(i,j);
        hs=hs_fIII(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fIII(i,j,1:3),gdxi_III(i,j,1:3));
        X(2)=dot(vt_fIII(i,j,1:3),gdeta_III(i,j,1:3));
        X(3)=ht_fIII(i,j);
        Y_fIII(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face IV
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fIV(i,j,1:3),gdxi_IV(i,j,1:3));
        ueta=dot(vt_fIV(i,j,1:3),gdeta_IV(i,j,1:3));
        G11=G11_fIV(i,j);
        G12=G12_fIV(i,j);
        G22=G22_fIV(i,j);
        h=ht_fIV(i,j);
        hs=hs_fIV(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fIV(i,j,1:3),gdxi_IV(i,j,1:3));
        X(2)=dot(vt_fIV(i,j,1:3),gdeta_IV(i,j,1:3));
        X(3)=ht_fIV(i,j);
        Y_fIV(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face V
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fV(i,j,1:3),gdxi_V(i,j,1:3));
        ueta=dot(vt_fV(i,j,1:3),gdeta_V(i,j,1:3));
        G11=G11_fV(i,j);
        G12=G12_fV(i,j);
        G22=G22_fV(i,j);
        h=ht_fV(i,j);
        hs=hs_fV(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fV(i,j,1:3),gdxi_V(i,j,1:3));
        X(2)=dot(vt_fV(i,j,1:3),gdeta_V(i,j,1:3));
        X(3)=ht_fV(i,j);
        Y_fV(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face VI
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fVI(i,j,1:3),gdxi_VI(i,j,1:3));
        ueta=dot(vt_fVI(i,j,1:3),gdeta_VI(i,j,1:3));
        G11=G11_fVI(i,j);
        G12=G12_fVI(i,j);
        G22=G22_fVI(i,j);
        h=ht_fVI(i,j);
        hs=hs_fVI(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        X(1)=dot(vt_fVI(i,j,1:3),gdxi_VI(i,j,1:3));
        X(2)=dot(vt_fVI(i,j,1:3),gdeta_VI(i,j,1:3));
        X(3)=ht_fVI(i,j);
        Y_fVI(i,j,1:3)=V_xi\(V_eta\X');
    end
end

%% *** Filter *************************************************************
[Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1)]=ftr72(Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1),n,nn);
[Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2)]=ftr72(Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2),n,nn);
[Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3)]=ftr72(Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3),n,nn);

%% *** back to classical variables ****************************************
% Face I
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fI(i,j,1:3),gdxi_I(i,j,1:3));
        ueta=dot(vt_fI(i,j,1:3),gdeta_I(i,j,1:3));
        G11=G11_fI(i,j);
        G12=G12_fI(i,j);
        G22=G22_fI(i,j);
        h=ht_fI(i,j);
        hs=hs_fI(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fI(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fI(i,j,1:3)=X(1).*gxi_I(i,j,1:3)+X(2)*geta_I(i,j,1:3);
        hf_fI(i,j)=X(3);
    end
end

% Face II
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fII(i,j,1:3),gdxi_II(i,j,1:3));
        ueta=dot(vt_fII(i,j,1:3),gdeta_II(i,j,1:3));
        G11=G11_fII(i,j);
        G12=G12_fII(i,j);
        G22=G22_fII(i,j);
        h=ht_fII(i,j);
        hs=hs_fII(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fII(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fII(i,j,1:3)=X(1).*gxi_II(i,j,1:3)+X(2)*geta_II(i,j,1:3);
        hf_fII(i,j)=X(3);
    end
end

% Face III
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fIII(i,j,1:3),gdxi_III(i,j,1:3));
        ueta=dot(vt_fIII(i,j,1:3),gdeta_III(i,j,1:3));
        G11=G11_fIII(i,j);
        G12=G12_fIII(i,j);
        G22=G22_fIII(i,j);
        h=ht_fIII(i,j);
        hs=hs_fIII(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fIII(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fIII(i,j,1:3)=X(1).*gxi_III(i,j,1:3)+X(2)*geta_III(i,j,1:3);
        hf_fIII(i,j)=X(3);
    end
end

% Face IV
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fIV(i,j,1:3),gdxi_IV(i,j,1:3));
        ueta=dot(vt_fIV(i,j,1:3),gdeta_IV(i,j,1:3));
        G11=G11_fIV(i,j);
        G12=G12_fIV(i,j);
        G22=G22_fIV(i,j);
        h=ht_fIV(i,j);
        hs=hs_fIV(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fIV(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fIV(i,j,1:3)=X(1).*gxi_IV(i,j,1:3)+X(2)*geta_IV(i,j,1:3);
        hf_fIV(i,j)=X(3);
    end
end

% Face V
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fV(i,j,1:3),gdxi_V(i,j,1:3));
        ueta=dot(vt_fV(i,j,1:3),gdeta_V(i,j,1:3));
        G11=G11_fV(i,j);
        G12=G12_fV(i,j);
        G22=G22_fV(i,j);
        h=ht_fV(i,j);
        hs=hs_fV(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fV(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fV(i,j,1:3)=X(1).*gxi_V(i,j,1:3)+X(2)*geta_V(i,j,1:3);
        hf_fV(i,j)=X(3);
    end
end

% Face VI
for i=1:nn
    for j=1:nn
        uxi=dot(vt_fVI(i,j,1:3),gdxi_VI(i,j,1:3));
        ueta=dot(vt_fVI(i,j,1:3),gdeta_VI(i,j,1:3));
        G11=G11_fVI(i,j);
        G12=G12_fVI(i,j);
        G22=G22_fVI(i,j);
        h=ht_fVI(i,j);
        hs=hs_fVI(i,j);
        
        % construction de V_xi
        V11=(G11^2*uxi^2 - G12*gp*h + G12*gp*hs + G12*G22*ueta^2 + G11*G12*ueta*uxi + G11*G22*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V21=-(G12^2*ueta^2 + G11^2*uxi^2 - G11*gp*h + G11*gp*hs + 2*G11*G12*ueta*uxi)/((h - hs)*(- uxi*G11^2 + uxi*G11*G12 - G22*ueta*G11 + ueta*G12^2));
        V31=1;
        V12=(G12*ueta + G11*uxi + (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V22=0;
        V32=1;
        V13=(G12*ueta + G11*uxi - (G11*gp*(h - hs))^(1/2))/(G11*(h - hs)) - (G12*ueta + G11*uxi)/(G11*(h - hs));
        V23=0;
        V33=1;
        
        V_xi=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        % construction de V_eta
        V11=-(G22^2*ueta^2 + G12^2*uxi^2 - G22*gp*h + G22*gp*hs + 2*G12*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V21=(G12*gp*hs - G12*gp*h + G12*G22*ueta^2 + G11*G12*uxi^2 + G12^2*ueta*uxi + G11*G22*ueta*uxi)/((uxi*G12^2 - G11*G22*uxi)*(h - hs));
        V31=1;
        V12=0;
        V22=(G22*ueta + G12*uxi + (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V32=1;
        V13=0;
        V23=(G22*ueta + G12*uxi - (G22*gp*(h - hs))^(1/2))/(G22*(h - hs)) - (G22*ueta + G12*uxi)/(G22*(h - hs));
        V33=1;
        
        V_eta=[V11 V12 V13;V21 V22 V23; V31 V32 V33];
        
        u=squeeze(Y_fVI(i,j,1:3));
        X(1:3)=V_eta*(V_xi*u);
        vf_fVI(i,j,1:3)=X(1).*gxi_VI(i,j,1:3)+X(2)*geta_VI(i,j,1:3);
        hf_fVI(i,j)=X(3);
    end
end
end