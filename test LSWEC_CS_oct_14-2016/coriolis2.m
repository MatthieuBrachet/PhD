function[flu_I,flu_II,flu_III,flu_IV,flu_V,flu_VI]=coriolis2(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global omega alpha
global gr_I gr_II gr_III gr_IV gr_V gr_VI;

%% coriolis

% PDT VECT
[n1,n2]=size(x_fI);
for i=1:n1
    for j=1:n2
        vect_fI(i,j,1:3)=cross(gr_I(i,j,1:3),vt_fI(i,j,1:3));
        vect_fII(i,j,1:3)=cross(gr_II(i,j,1:3),vt_fII(i,j,1:3));
        vect_fIII(i,j,1:3)=cross(gr_III(i,j,1:3),vt_fIII(i,j,1:3));
        vect_fIV(i,j,1:3)=cross(gr_IV(i,j,1:3),vt_fIV(i,j,1:3));
        vect_fV(i,j,1:3)=cross(gr_V(i,j,1:3),vt_fV(i,j,1:3));
        vect_fVI(i,j,1:3)=cross(gr_VI(i,j,1:3),vt_fVI(i,j,1:3));
    end
end

%% ASSEMBLAGE
[lambda_fI, teta_fI, aaa]=cart2sph(x_fI, y_fI, z_fI);
[lambda_fII, teta_fII, aaa]=cart2sph(x_fII, y_fII, z_fII);
[lambda_fIII, teta_fIII, aaa]=cart2sph(x_fIII, y_fIII, z_fIII);
[lambda_fIV, teta_fIV, aaa]=cart2sph(x_fIV, y_fIV, z_fIV);
[lambda_fV, teta_fV, aaa]=cart2sph(x_fV, y_fV, z_fV);
[lambda_fVI, teta_fVI, aaa]=cart2sph(x_fVI, y_fVI, z_fVI);

for kk=1:3
    f_I=2.*omega.*sin(teta_fI);
    flu_I(:,:,kk)=f_I.*vect_fI(:,:,kk);
    
    f_II=2.*omega.*sin(teta_fII);
    flu_II(:,:,kk)=f_II.*vect_fII(:,:,kk);
    
    f_III=2.*omega.*sin(teta_fIII);
    flu_III(:,:,kk)=f_III.*vect_fIII(:,:,kk);
    
    f_IV=2.*omega.*sin(teta_fIV);
    flu_IV(:,:,kk)=f_IV.*vect_fIV(:,:,kk);
    
    f_V=2.*omega.*sin(teta_fV);
    flu_V(:,:,kk)=f_V.*vect_fV(:,:,kk);
    
    f_VI=2.*omega.*sin(teta_fVI);
    flu_VI(:,:,kk)=f_VI.*vect_fVI(:,:,kk);
end

