clc; clear all; close all;
%% test sur FUN5

clc; clear all; close all; format short;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp u0 radius omega

test=1;
video = 'no';
save = 0;
opt_ftr=10;
n=20;
mod72

x=linspace(-radius,radius,100);
% y=fun5(x);
% yb=fun5b(x);

y=fun7(x);
yb=fun7b(x);
