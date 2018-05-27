clc; clear all; close all;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=2;
make_cs_grid(N);

% Computation of the weights with the basic formula
weights=dxi*deta*dga;

 
 
nhs1=5;mhs1=0;

for nhs2=0:nhs1
    if nhs2<nhs1
        mhs2_max=nhs2;
    elseif nhs2==nhs1
        mhs2_max=mhs1;
    else
        error('erreur');
    end
    
    for mhs2=0:mhs2_max

        fun1fI=conj(sph(nhs1,mhs1,x_fI,y_fI,z_fI));
        fun1fII=conj(sph(nhs1,mhs1,x_fII,y_fII,z_fII));
        fun1fIII=conj(sph(nhs1,mhs1,x_fIII,y_fIII,z_fIII));
        fun1fIV=conj(sph(nhs1,mhs1,x_fIV,y_fIV,z_fIV));
        fun1fV=conj(sph(nhs1,mhs1,x_fV,y_fV,z_fV));
        fun1fVI=conj(sph(nhs1,mhs1,x_fVI,y_fVI,z_fVI));

        fun2fI=sph(nhs2,mhs2,x_fI,y_fI,z_fI);
        fun2fII=sph(nhs2,mhs2,x_fII,y_fII,z_fII);
        fun2fIII=sph(nhs2,mhs2,x_fIII,y_fIII,z_fIII);
        fun2fIV=sph(nhs2,mhs2,x_fIV,y_fIV,z_fIV);
        fun2fV=sph(nhs2,mhs2,x_fV,y_fV,z_fV);
        fun2fVI=sph(nhs2,mhs2,x_fVI,y_fVI,z_fVI);

        % APPROXIMATE QSCALAR PRODUCT ON THE SPHERE
        [scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
            intsca_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
            fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI);
        
        disp(['n''       = ' num2str(nhs2)])
        disp(['m''       = ' num2str(mhs2)])
        disp(['scal     = ' num2str(scaf1f2)])
        disp(['abs(scal) = ' num2str(abs(scaf1f2))])
        disp('**********************')
    end
end