clc;
clear all;
% fun2 is spherical harmonic arounf (Oz),
% fun2 is spherical harmonic arounf (Oy).

global n nn
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=32;
make_cs_grid(N);

% Computation of the weight with the basic formula
weights=dxi*deta*dga;


nhs1=23; mhs1=12;
nhs2=22; mhs2=12;

fun1fI=conj(sph(nhs1,mhs1,x_fI,y_fI,z_fI));
fun1fII=conj(sph(nhs1,mhs1,x_fII,y_fII,z_fII));
fun1fIII=conj(sph(nhs1,mhs1,x_fIII,y_fIII,z_fIII));
fun1fIV=conj(sph(nhs1,mhs1,x_fIV,y_fIV,z_fIV));
fun1fV=conj(sph(nhs1,mhs1,x_fV,y_fV,z_fV));
fun1fVI=conj(sph(nhs1,mhs1,x_fVI,y_fVI,z_fVI));

fun2fI=sph(nhs2,mhs2,x_fI,-z_fI,y_fI);
fun2fII=sph(nhs2,mhs2,x_fII,-z_fII,y_fII);
fun2fIII=sph(nhs2,mhs2,x_fIII,-z_fIII,y_fIII);
fun2fIV=sph(nhs2,mhs2,x_fIV,-z_fIV,y_fIV);
fun2fV=sph(nhs2,mhs2,x_fV,-z_fV,y_fV);
fun2fVI=sph(nhs2,mhs2,x_fVI,-z_fVI,y_fVI);

figure(1)
plot_cs11(n,nn,real(fun1fI),real(fun1fII),real(fun1fIII),real(fun1fIV),real(fun1fV),real(fun1fVI))

figure(2)
plot_cs11(n,nn,real(fun2fI),real(fun2fII),real(fun2fIII),real(fun2fIV),real(fun2fV),real(fun2fVI))

% APPROXIMATE QSCALAR PRODUCT ON THE SPHERE
[scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
    intsca_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
    fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI);

scaf1f2

fig_placier