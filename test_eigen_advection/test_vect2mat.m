clc; clear all; close all;
global n nn
N=32;
n=N-1;
nn=N+1;
Nb=nn*nn;


%% données
funfI=rand(nn,nn);
funfII=rand(nn,nn);
funfIII=rand(nn,nn);
funfIV=rand(nn,nn);
funfV=rand(nn,nn);
funfVI=rand(nn,nn);

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ds101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

funfI=floor(funfI);
funfII=floor(funfII);
funfIII=floor(funfIII);
funfIV=floor(funfIV);
funfV=floor(funfV);
funfVI=floor(funfVI);

%% fun to vect.
[ fun ] = fun2vect( funfI,funfII,funfIII,funfIV,funfV,funfVI );
[ funI, funII, funIII, funIV, funV, funVI ] = vect2fun( fun );

eI=funI-funfI;
eII=funII-funfII;
eIII=funIII-funfIII;
eIV=funIV-funfIV;
eV=funV-funfV;
eVI=funVI-funfVI;

str='infty';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(eI,eII,eIII,eIV,eV,eVI,n,nn,str);
nrmg

Nl=length(fun);
Nl-(6*N^2+2)