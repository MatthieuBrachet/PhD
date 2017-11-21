clc; clear all; close all; format shorte
global n nn
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global scheme

scheme='compact4';
n=63;
%mod_bis;
mod101

[mfunfI,gr_fI] = fun2(x_fI,y_fI,z_fI);
[mfunfII,gr_fII] = fun2(x_fII,y_fII,z_fII);
[mfunfIII,gr_fIII] = fun2(x_fIII,y_fIII,z_fIII);
[mfunfIV,gr_fIV] = fun2(x_fIV,y_fIV,z_fIV);
[mfunfV,gr_fV] = fun2(x_fV,y_fV,z_fV);
[mfunfVI,gr_fVI] = fun2(x_fVI,y_fVI,z_fVI);

MM=max(max(max(abs([gr_fI gr_fII gr_fIII gr_fIV gr_fV gr_fVI]))));

%[grad_fI,grad_fII,grad_fIII,grad_fIV,grad_fV,grad_fVI]=grad1000(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);
[grad_fI,grad_fII,grad_fIII,grad_fIV,grad_fV,grad_fVI]=gr101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);
err=max(max(max(abs([grad_fI-gr_fI; grad_fII-gr_fII; grad_fIII-gr_fIII; grad_fIV-gr_fIV; grad_fV-gr_fV; grad_fVI-gr_fVI]))))./MM
% 
% err_I=((grad_fI(:,:,1)-gr_fI(:,:,1)).^2+(grad_fI(:,:,2)-gr_fI(:,:,2)).^2+(grad_fI(:,:,3)-gr_fI(:,:,3)).^2)./MM.^2;
% err_II=((grad_fII(:,:,1)-gr_fII(:,:,1)).^2+(grad_fII(:,:,2)-gr_fII(:,:,2)).^2+(grad_fII(:,:,3)-gr_fII(:,:,3)).^2)./MM.^2;
% err_III=((grad_fIII(:,:,1)-gr_fIII(:,:,1)).^2+(grad_fIII(:,:,2)-gr_fIII(:,:,2)).^2+(grad_fIII(:,:,3)-gr_fIII(:,:,3)).^2)./MM.^2;
% err_IV=((grad_fIV(:,:,1)-gr_fIV(:,:,1)).^2+(grad_fIV(:,:,2)-gr_fIV(:,:,2)).^2+(grad_fIV(:,:,3)-gr_fIV(:,:,3)).^2)./MM.^2;
% err_V=((grad_fV(:,:,1)-gr_fV(:,:,1)).^2+(grad_fV(:,:,2)-gr_fV(:,:,2)).^2+(grad_fV(:,:,3)-gr_fV(:,:,3)).^2)./MM.^2;
% err_VI=((grad_fVI(:,:,1)-gr_fVI(:,:,1)).^2+(grad_fVI(:,:,2)-gr_fVI(:,:,2)).^2+(grad_fVI(:,:,3)-gr_fVI(:,:,3)).^2)./MM.^2;
% 
% 
% figure(1)
% plot_cs100(n,nn,err_I,err_II,err_III,err_IV,err_V,err_VI)
% title('norm(error)')
% colorbar
% 
% figure(2)
% plot_cs15(n,nn,mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI)
% colorbar
% title('function')
% 
% e_I=max(abs(grad_fI-gr_fI),[],3)./MM;
% e_II=max(abs(grad_fII-gr_fII),[],3)./MM;
% e_III=max(abs(grad_fIII-gr_fIII),[],3)./MM;
% e_IV=max(abs(grad_fIV-gr_fIV),[],3)./MM;
% e_V=max(abs(grad_fV-gr_fV),[],3)./MM;
% e_VI=max(abs(grad_fVI-gr_fVI),[],3)./MM;
% 
% hFig=figure(3)
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [50 50 1000 500])
% plot_cs100(n,nn,e_I,e_II,e_III,e_IV,e_V,e_VI)
% colorbar
% 
% 
% fig_placier