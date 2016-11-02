function [ hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI, vf_fI, vf_fII, vf_fIII, vf_fIV, vf_fV, vf_fVI ]...
    = ftrcar72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
%% *** caracteristic filter
% caracteristic filter for the LSWE problem : 
%
%   dv/dt + f * k vect v + gp*grad(h) = 0
%   dh/dt + H * div(v) = 0
%
% created the 10/10/2016 by Matthieu Brachet.

global gp hp
global nn n
% -----------------------------------
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global gdxi_I gdxi_II gdxi_III gdxi_IV gdxi_V gdxi_VI;
global gdeta_I gdeta_II gdeta_III gdeta_IV gdeta_V gdeta_VI;
% -----------------------------------

%% *** caracteristics variables *******************************************
% Face I
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        
        X(1)=dot(vt_fI(i,j,1:3),gdxi_I(i,j,1:3));
        X(2)=dot(vt_fI(i,j,1:3),gdeta_I(i,j,1:3));
        X(3)=ht_fI(i,j);
        Y_fI(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face II
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        
        X(1)=dot(vt_fII(i,j,1:3),gdxi_II(i,j,1:3));
        X(2)=dot(vt_fII(i,j,1:3),gdeta_II(i,j,1:3));
        X(3)=ht_fII(i,j);
        Y_fII(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face III
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        
        X(1)=dot(vt_fIII(i,j,1:3),gdxi_III(i,j,1:3));
        X(2)=dot(vt_fIII(i,j,1:3),gdeta_III(i,j,1:3));
        X(3)=ht_fIII(i,j);
        Y_fIII(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face IV
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        
        X(1)=dot(vt_fIV(i,j,1:3),gdxi_IV(i,j,1:3));
        X(2)=dot(vt_fIV(i,j,1:3),gdeta_IV(i,j,1:3));
        X(3)=ht_fIV(i,j);
        Y_fIV(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face V
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        
        X(1)=dot(vt_fV(i,j,1:3),gdxi_V(i,j,1:3));
        X(2)=dot(vt_fV(i,j,1:3),gdeta_V(i,j,1:3));
        X(3)=ht_fV(i,j);
        Y_fV(i,j,1:3)=V_xi\(V_eta\X');
    end
end

% Face VI
for i=1:nn
    for j=1:nn
        h=hp;
        
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
        h=hp;
        
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
        h=hp;
        
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
        h=hp;
        
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
        h=hp;
        
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
        h=hp;
        
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
        h=hp;
        
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