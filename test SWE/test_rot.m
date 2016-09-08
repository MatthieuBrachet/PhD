%% test rotational
clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr
global teta0 teta1


teta0=-3*pi/16;
teta1=3*pi/16;
opt_ftr=0;
n=100;
mod72;

[vt_fI,rotvt_fI] = sol_exacte_rot(x_fI,y_fI,z_fI);
[vt_fII,rotvt_fII] = sol_exacte_rot(x_fII,y_fII,z_fII);
[vt_fIII,rotvt_fIII] = sol_exacte_rot(x_fIII,y_fIII,z_fIII);
[vt_fIV,rotvt_fIV] = sol_exacte_rot(x_fIV,y_fIV,z_fIV);
[vt_fV,rotvt_fV] = sol_exacte_rot(x_fV,y_fV,z_fV);
[vt_fVI,rotvt_fVI] = sol_exacte_rot(x_fVI,y_fVI,z_fVI);

[rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

err_I=abs(rotvt_fI-rot_fI);
err_II=abs(rotvt_fII-rot_fII);
err_III=abs(rotvt_fIII-rot_fIII);
err_IV=abs(rotvt_fIV-rot_fIV);
err_V=abs(rotvt_fV-rot_fV);
err_VI=abs(rotvt_fVI-rot_fVI);

e_I=max(max(max(err_I)))
e_II=max(max(max(err_II)))
e_III=max(max(max(err_III)))
e_IV=max(max(max(err_IV)))
e_V=max(max(max(err_V)))
e_VI=max(max(max(err_VI)))

err=max([e_I e_II e_III e_IV e_V e_VI])

% for p=1:3
%     m_I(p)=mean(mean(rotvt_fI(:,:,p)./rot_fI(:,:,p)));
%     m_II(p)=mean(mean(rotvt_fII(:,:,p)./rot_fII(:,:,p)));
%     m_III(p)=mean(mean(rotvt_fIII(:,:,p)./rot_fIII(:,:,p)));
%     m_IV(p)=mean(mean(rotvt_fIV(:,:,p)./rot_fIV(:,:,p)));
%     m_V(p)=mean(mean(rotvt_fV(:,:,p)./rot_fV(:,:,p)));
%     m_VI(p)=mean(mean(rotvt_fVI(:,:,p)./rot_fVI(:,:,p)));
% end
% 
% m_I
% m_II
% m_III
% m_IV
% m_V
% m_VI


% for p=1:3
%     figure(p)
%     plot_cs11(n,nn,rotvt_fI(:,:,p)./rot_fI(:,:,p),rotvt_fII(:,:,p)./rot_fII(:,:,p),rotvt_fIII(:,:,p)./rot_fIII(:,:,p),rotvt_fIV(:,:,p)./rot_fIV(:,:,p),rotvt_fV(:,:,p)./rot_fV(:,:,p),rotvt_fVI(:,:,p)./rot_fVI(:,:,p))
%     title('quotient exact/approx')
%     colorbar
% end
% 
% for p=1:3
%     figure(p+3)
%     subplot(1,3,1)
%     plot_cs11(n,nn,rotvt_fI(:,:,p),rotvt_fII(:,:,p),rotvt_fIII(:,:,p),rotvt_fIV(:,:,p),rotvt_fV(:,:,p),rotvt_fVI(:,:,p))
%     title('exact sol.')
%     colorbar
% 
%     subplot(1,3,2)
%     plot_cs11(n,nn,rot_fI(:,:,p),rot_fII(:,:,p),rot_fIII(:,:,p),rot_fIV(:,:,p),rot_fV(:,:,p),rot_fVI(:,:,p))
%     title('approx sol.')
%     colorbar
%     
%     subplot(1,3,3)
%     plot_cs11(n,nn,err_I(:,:,p),err_II(:,:,p),err_III(:,:,p),err_IV(:,:,p),err_V(:,:,p),err_VI(:,:,p))
%     title('ERROR')
%     colorbar
% end


