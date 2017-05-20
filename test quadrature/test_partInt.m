clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme filtre test
global radius dga;

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
test=1;

n=127; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
mod101
disp('mod74 : ok')


funfI=ones(size(x_fI));
funfII=zeros(size(x_fI));
funfIII=zeros(size(x_fI));
funfIV=zeros(size(x_fI));
funfV=zeros(size(x_fI));
funfVI=zeros(size(x_fI));

area=4*pi*radius.^2;

str='int'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);




figure(1)
plot_cs11(n,nn,dga/(radius.^2),dga/(radius.^2),dga/(radius.^2),dga/(radius.^2),dga/(radius.^2),dga/(radius.^2))
