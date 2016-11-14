function [ hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI, vf_fI, vf_fII, vf_fIII, vf_fIV, vf_fV, vf_fVI ]...
    = ftrcar74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
%% *** caracteristic filter
% caracteristic filter for the SWE problem : 
%
%   dv/dt + (v*nabla)*v + f * k vect v + gp*grad(h) = 0
%   dh/dt + div(h*v) = 0
%
% created the 26/10/2016 by Matthieu Brachet.

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

% face I
[n1,n2]=size(ht_fI);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fI(i,j,1:3),gxi_I(i,j,1:3));
        X(2)=dot(vt_fI(i,j,1:3),geta_I(i,j,1:3));
        X(3)=ht_fI(i,j);
        
        h=ht_fI(i,j);
        u=gxi_I(i,j,1:3);
        v=geta_I(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fI(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end

% face II
[n1,n2]=size(ht_fII);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fII(i,j,1:3),gxi_II(i,j,1:3));
        X(2)=dot(vt_fII(i,j,1:3),geta_II(i,j,1:3));
        X(3)=ht_fII(i,j);
        
        h=ht_fII(i,j);
        u=gxi_II(i,j,1:3);
        v=geta_II(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fII(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end

% face III
[n1,n2]=size(ht_fIII);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fIII(i,j,1:3),gxi_III(i,j,1:3));
        X(2)=dot(vt_fIII(i,j,1:3),geta_III(i,j,1:3));
        X(3)=ht_fIII(i,j);
        
        h=ht_fIII(i,j);
        u=gxi_III(i,j,1:3);
        v=geta_III(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fIII(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end

% face IV
[n1,n2]=size(ht_fIV);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fIV(i,j,1:3),gxi_IV(i,j,1:3));
        X(2)=dot(vt_fIV(i,j,1:3),geta_IV(i,j,1:3));
        X(3)=ht_fIV(i,j);
        
        h=ht_fIV(i,j);
        u=gxi_IV(i,j,1:3);
        v=geta_IV(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fIV(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end

% face V
[n1,n2]=size(ht_fV);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fV(i,j,1:3),gxi_V(i,j,1:3));
        X(2)=dot(vt_fV(i,j,1:3),geta_V(i,j,1:3));
        X(3)=ht_fV(i,j);
        
        h=ht_fV(i,j);
        u=gxi_V(i,j,1:3);
        v=geta_V(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fV(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end

% face VI
[n1,n2]=size(ht_fVI);
for i=1:n1
    for j=1:n2
        X(1)=dot(vt_fVI(i,j,1:3),gxi_VI(i,j,1:3));
        X(2)=dot(vt_fVI(i,j,1:3),geta_VI(i,j,1:3));
        X(3)=ht_fVI(i,j);
        
        h=ht_fVI(i,j);
        u=gxi_VI(i,j,1:3);
        v=geta_VI(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );

        Y_fVI(i,j,1:3)=V_xi\(V_eta\transpose(X));
    end
end


%% *** Filter *************************************************************
figure(100)
plot_cs11(n,nn,Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1));

[Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1)]=...
    ftr74(Y_fI(:,:,1),Y_fII(:,:,1),Y_fIII(:,:,1),Y_fIV(:,:,1),Y_fV(:,:,1),Y_fVI(:,:,1),n,nn);

figure(101)
plot_cs11(n,nn,Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2));

[Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2)]=...
    ftr74(Y_fI(:,:,2),Y_fII(:,:,2),Y_fIII(:,:,2),Y_fIV(:,:,2),Y_fV(:,:,2),Y_fVI(:,:,2),n,nn);

figure(102)
plot_cs11(n,nn,Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3));

[Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3)]=...
    ftr74(Y_fI(:,:,3),Y_fII(:,:,3),Y_fIII(:,:,3),Y_fIV(:,:,3),Y_fV(:,:,3),Y_fVI(:,:,3),n,nn);

%% *** back to classical variables ****************************************

% Face I
[n1,n2]=size(ht_fI);
for i=1:n1
    for j=1:n2
        h=ht_fI(i,j);
        u=gxi_I(i,j,1:3);
        v=geta_I(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fI(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fI(i,j,1:3)=uxi*gdxi_I(i,j,1:3)+ueta*gdeta_I(i,j,1:3);
        hf_fI(i,j)=Y(3);
    end
end

% Face II
[n1,n2]=size(ht_fII);
for i=1:n1
    for j=1:n2
        h=ht_fII(i,j);
        u=gxi_II(i,j,1:3);
        v=geta_II(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fII(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fII(i,j,1:3)=uxi*gdxi_II(i,j,1:3)+ueta*gdeta_II(i,j,1:3);
        hf_fII(i,j)=Y(3);
    end
end

% Face III
[n1,n2]=size(ht_fIII);
for i=1:n1
    for j=1:n2
        h=ht_fIII(i,j);
        u=gxi_III(i,j,1:3);
        v=geta_III(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fIII(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fIII(i,j,1:3)=uxi*gdxi_III(i,j,1:3)+ueta*gdeta_III(i,j,1:3);
        hf_fIII(i,j)=Y(3);
    end
end

% Face IV
[n1,n2]=size(ht_fIV);
for i=1:n1
    for j=1:n2
        h=ht_fIV(i,j);
        u=gxi_IV(i,j,1:3);
        v=geta_IV(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fIV(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fIV(i,j,1:3)=uxi*gdxi_IV(i,j,1:3)+ueta*gdeta_IV(i,j,1:3);
        hf_fIV(i,j)=Y(3);
    end
end

% Face V
[n1,n2]=size(ht_fV);
for i=1:n1
    for j=1:n2
        h=ht_fV(i,j);
        u=gxi_V(i,j,1:3);
        v=geta_V(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fV(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fV(i,j,1:3)=uxi*gdxi_V(i,j,1:3)+ueta*gdeta_V(i,j,1:3);
        hf_fV(i,j)=Y(3);
    end
end

% Face VI
[n1,n2]=size(ht_fVI);
for i=1:n1
    for j=1:n2
        h=ht_fVI(i,j);
        u=gxi_VI(i,j,1:3);
        v=geta_VI(i,j,1:3);
        [ V_xi, V_eta ] = riemann74( h, u, v );
        
        X=squeeze(Y_fVI(i,j,1:3));
        Y=V_eta*(V_xi*X);
        
        uxi=Y(1);
        ueta=Y(2);
        vf_fVI(i,j,1:3)=uxi*gdxi_VI(i,j,1:3)+ueta*gdeta_VI(i,j,1:3);
        hf_fVI(i,j)=Y(3);
    end
end
