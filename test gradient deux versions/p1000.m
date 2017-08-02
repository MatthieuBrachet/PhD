clc; clear all; close all;
global n nn
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;


n=31;
mod_bis;


[mfunfI] = fun(x_fI,y_fI,z_fI);
[mfunfII] = fun(x_fII,y_fII,z_fII);
[mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
[mfunfIV] = fun(x_fIV,y_fIV,z_fIV);
[mfunfV] = fun(x_fV,y_fV,z_fV);
[mfunfVI] = fun(x_fVI,y_fVI,z_fVI);

[grad_fI,grad_fII,grad_fIII,grad_fIV,grad_fV,grad_fVI]=grad1000(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

mod101;
[grad2_fI,grad2_fII,grad2_fIII,grad2_fIV,grad2_fV,grad2_fVI]=gr101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

for ppp=1:3
    err_I=grad_fI(:,:,ppp)-grad2_fI(:,:,ppp);
    err_II=grad_fII(:,:,ppp)-grad2_fII(:,:,ppp);
    err_III=grad_fIII(:,:,ppp)-grad2_fIII(:,:,ppp);
    err_IV=grad_fIV(:,:,ppp)-grad2_fIV(:,:,ppp);
    err_V=grad_fV(:,:,ppp)-grad2_fV(:,:,ppp);
    err_VI=grad_fVI(:,:,ppp)-grad2_fVI(:,:,ppp);
    
    figure(ppp)
    plot_cs100(n,nn,err_I,err_II,err_III,err_IV,err_V,err_VI)
    colorbar
end

fig_placier

