clc; clear all; close all;
format long

global n nn radius
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme filtre

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

n=1023; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
mod101
disp('mod74 : ok')

nhs=10;
mhs=4;
fun_I=sph(nhs,mhs,x_fI,y_fI,z_fI);
fun_II=sph(nhs,mhs,x_fII,y_fII,z_fII);
fun_III=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
fun_IV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
fun_V=sph(nhs,mhs,x_fV,y_fV,z_fV);
fun_VI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);

str='int'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg1)

str='test'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg2)

str='simpson'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg3]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg3)




