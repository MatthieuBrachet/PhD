clc; clear all; close all;
format long;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
N=1;
n=N-1;nn=N+1;
plot_cs1_mesh(n,nn);
%%%%%%%%%%%%%%%%%%%%
figure(2);
N=2;
n=N-1;nn=N+1;
plot_cs1_mesh(n,nn);
%%%%%%%%%%%%%%%%%%%%
figure(3);
N=4;
n=N-1;nn=N+1;
plot_cs1_mesh(n,nn);
%%%%%%%%%%%%%%%%%%%%%%
figure(4);
N=6;
n=N-1,nn=N+1;
plot_cs1_mesh(n,nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=2;
make_cs_grid(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the weights with the basic formula
% weights=dxi*deta*dga;
nhs1_a=50;
nhs1_b=50
%mhs1=0;
% mhs1=nhs1;
for nhsk=nhs1_a:nhs1_b,
    % for mhs2=0:mhs1,
    for mhsk=0:nhsk,
        fun_fI=sph(nhsk,mhsk,x_fI,y_fI,z_fI);
        fun_fII=sph(nhsk,mhsk,x_fII,y_fII,z_fII);
        fun_fIII=sph(nhsk,mhsk,x_fIII,y_fIII,z_fIII);
        fun_fIV=sph(nhsk,mhsk,x_fIV,y_fIV,z_fIV);
        fun_fV=sph(nhsk,mhsk,x_fV,y_fV,z_fV);
        fun_fVI=sph(nhsk,mhsk,x_fVI,y_fVI,z_fVI);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['n''       = ' num2str(nhsk)])
        disp(['m''       = ' num2str(mhsk)])
        disp(['max_fI   = ' num2str(max(max(abs(fun_fI))))])
        disp(['max_fII  = ' num2str(max(max(abs(fun_fII))))])
        disp(['max_fIII = ' num2str(max(max(abs(fun_fIII))))])
        disp(['max_fIV  = ' num2str(max(max(abs(fun_fIV))))])
        disp(['max_fV   = ' num2str(max(max(abs(fun_fV))))])
        disp(['max_fVI  = ' num2str(max(max(abs(fun_fVI))))])
        disp('**********************')
    end
end
% break
% for nhs2=0:nhs1
%     if nhs2<nhs1
%         mhs2_max=nhs2;
%     elseif nhs2==nhs1
%         mhs2_max=mhs1;
%     else
%         error('erreur');
%     end
%     
%     for mhs2=0:mhs2_max
% 
%         fun1fI=conj(sph(nhs1,mhs1,x_fI,y_fI,z_fI));
%         fun1fII=conj(sph(nhs1,mhs1,x_fII,y_fII,z_fII));
%         fun1fIII=conj(sph(nhs1,mhs1,x_fIII,y_fIII,z_fIII));
%         fun1fIV=conj(sph(nhs1,mhs1,x_fIV,y_fIV,z_fIV));
%         fun1fV=conj(sph(nhs1,mhs1,x_fV,y_fV,z_fV));
%         fun1fVI=conj(sph(nhs1,mhs1,x_fVI,y_fVI,z_fVI));
% 
%         fun2fI=sph(nhs2,mhs2,x_fI,y_fI,z_fI);
%         fun2fII=sph(nhs2,mhs2,x_fII,y_fII,z_fII);
%         fun2fIII=sph(nhs2,mhs2,x_fIII,y_fIII,z_fIII);
%         fun2fIV=sph(nhs2,mhs2,x_fIV,y_fIV,z_fIV);
%         fun2fV=sph(nhs2,mhs2,x_fV,y_fV,z_fV);
%         fun2fVI=sph(nhs2,mhs2,x_fVI,y_fVI,z_fVI);
% 
%         % APPROXIMATE QSCALAR PRODUCT ON THE SPHERE
%         [scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
%             intsca_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
%             fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI);
%         
%         disp(['n''       = ' num2str(nhs2)])
%         disp(['m''       = ' num2str(mhs2)])
%         disp(['scal     = ' num2str(scaf1f2)])
%         disp(['abs(scal) = ' num2str(abs(scaf1f2))])
%         disp('**********************')
%     end
% end

fig_placier