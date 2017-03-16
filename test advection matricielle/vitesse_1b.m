function [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(t)

global radius tmax;
global alphad u0;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

global coef

% test de Nair et Machenhauer
global teta_p lambda_p
global gamma rho0

% test de Nair et Jablonowski
global teta0 lambda0

% test de Nair et Lauritzen (deformational test case)
global kk

a=radius;

if coef == 0

%% ------------------------------------------------------------------------
%%                       CALCUL DE LA VITESSE
%% ------------------------------------------------------------------------

%% solid body - test de williamson
[lambda,teta,radius1]=cart2sph(x_fI,y_fI,z_fI);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_I=ulambda.*elambda_x + vteta.*eteta_x;
vity_I=ulambda.*elambda_y + vteta.*eteta_y;
vitz_I=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE II---------------------
[lambda,teta,radius1]=cart2sph(x_fII,y_fII,z_fII);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_II=ulambda.*elambda_x + vteta.*eteta_x;
vity_II=ulambda.*elambda_y + vteta.*eteta_y;
vitz_II=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE III---------------------
[lambda,teta,radius1]=cart2sph(x_fIII,y_fIII,z_fIII);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_III=ulambda.*elambda_x + vteta.*eteta_x;
vity_III=ulambda.*elambda_y + vteta.*eteta_y;
vitz_III=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE IV---------------------
[lambda,teta,radius1]=cart2sph(x_fIV,y_fIV,z_fIV);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_IV=ulambda.*elambda_x + vteta.*eteta_x;
vity_IV=ulambda.*elambda_y + vteta.*eteta_y;
vitz_IV=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE V---------------------
[lambda,teta,radius1]=cart2sph(x_fV,y_fV,z_fV);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_V=ulambda.*elambda_x + vteta.*eteta_x;
vity_V=ulambda.*elambda_y + vteta.*eteta_y;
vitz_V=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE VI---------------------
[lambda,teta,radius1]=cart2sph(x_fVI,y_fVI,z_fVI);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_VI=ulambda.*elambda_x + vteta.*eteta_x;
vity_VI=ulambda.*elambda_y + vteta.*eteta_y;
vitz_VI=ulambda.*elambda_z + vteta.*eteta_z;

elseif coef == 1
%% Nair and Machenhauer
% FACE I

    [lambda, teta, ~]=cart2sph(x_fI,y_fI,z_fI);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_I=ulambda.*elambda_x + vteta.*eteta_x;
    vity_I=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_I=ulambda.*elambda_z + vteta.*eteta_z;

% FACE II

    [lambda, teta, ~]=cart2sph(x_fII,y_fII,z_fII);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_II=ulambda.*elambda_x + vteta.*eteta_x;
    vity_II=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_II=ulambda.*elambda_z + vteta.*eteta_z;

% FACE III

    [lambda, teta, ~]=cart2sph(x_fIII,y_fIII,z_fIII);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_III=ulambda.*elambda_x + vteta.*eteta_x;
    vity_III=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_III=ulambda.*elambda_z + vteta.*eteta_z;

% FACE IV

    [lambda, teta, ~]=cart2sph(x_fIV,y_fIV,z_fIV);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_IV=ulambda.*elambda_x + vteta.*eteta_x;
    vity_IV=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_IV=ulambda.*elambda_z + vteta.*eteta_z;

% FACE V

    [lambda, teta, ~]=cart2sph(x_fV,y_fV,z_fV);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_V=ulambda.*elambda_x + vteta.*eteta_x;
    vity_V=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_V=ulambda.*elambda_z + vteta.*eteta_z;

% FACE VI

    [lambda, teta, ~]=cart2sph(x_fVI,y_fVI,z_fVI);
    [ ~, teta_prime ] = rotated_coord( lambda, teta );
    rho=rho0*cos(teta_prime);
    V=u0*1.5*sqrt(3)*sech(rho).^2.*tanh(rho);
    logi=(rho~=0);
    w_r=V./(radius*rho).*logi;
    
    ulambda=a*w_r.*(sin(teta_p).*cos(teta)-cos(teta_p).*cos(lambda-lambda_p).*sin(teta));
    vteta=a*w_r.*(cos(teta_p).*sin(lambda-lambda_p));
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);
    
    % en coordonnées cartesiennes
    vitx_VI=ulambda.*elambda_x + vteta.*eteta_x;
    vity_VI=ulambda.*elambda_y + vteta.*eteta_y;
    vitz_VI=ulambda.*elambda_z + vteta.*eteta_z;

elseif coef == 2
%% Nair and Jablonowski
    [vitx_I   , vity_I   , vitz_I  ] = nair_jablonowski(x_fI  ,y_fI   ,z_fI  ,t);
    [vitx_II  , vity_II  , vitz_II ] = nair_jablonowski(x_fII ,y_fII  ,z_fII ,t);
    [vitx_III , vity_III , vitz_III] = nair_jablonowski(x_fIII,y_fIII ,z_fIII,t);
    [vitx_IV  , vity_IV  , vitz_IV ] = nair_jablonowski(x_fIV ,y_fIV  ,z_fIV ,t);
    [vitx_V   , vity_V   , vitz_V  ] = nair_jablonowski(x_fV  ,y_fV   ,z_fV  ,t);
    [vitx_VI  , vity_VI  , vitz_VI ] = nair_jablonowski(x_fVI ,y_fVI  ,z_fVI ,t);
    
elseif coef == 3
    
%% slotted cylinder - test de zalesak

% --------------- FACE I---------------------
[lambda,teta,radius1]=cart2sph(x_fI,y_fI,z_fI);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_I=ulambda.*elambda_x + vteta.*eteta_x;
vity_I=ulambda.*elambda_y + vteta.*eteta_y;
vitz_I=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE II---------------------
[lambda,teta,radius1]=cart2sph(x_fII,y_fII,z_fII);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_II=ulambda.*elambda_x + vteta.*eteta_x;
vity_II=ulambda.*elambda_y + vteta.*eteta_y;
vitz_II=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE III---------------------
[lambda,teta,radius1]=cart2sph(x_fIII,y_fIII,z_fIII);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_III=ulambda.*elambda_x + vteta.*eteta_x;
vity_III=ulambda.*elambda_y + vteta.*eteta_y;
vitz_III=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE IV---------------------
[lambda,teta,radius1]=cart2sph(x_fIV,y_fIV,z_fIV);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_IV=ulambda.*elambda_x + vteta.*eteta_x;
vity_IV=ulambda.*elambda_y + vteta.*eteta_y;
vitz_IV=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE V---------------------
[lambda,teta,radius1]=cart2sph(x_fV,y_fV,z_fV);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_V=ulambda.*elambda_x + vteta.*eteta_x;
vity_V=ulambda.*elambda_y + vteta.*eteta_y;
vitz_V=ulambda.*elambda_z + vteta.*eteta_z;
% --------------- FACE VI---------------------
[lambda,teta,radius1]=cart2sph(x_fVI,y_fVI,z_fVI);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z=0;
eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
%
ulambda=u0*(cos(alphad)*cos(teta)+sin(alphad)*sin(teta).*cos(lambda));
vteta= -u0*sin(alphad)*sin(lambda);
%
vitx_VI=ulambda.*elambda_x + vteta.*eteta_x;
vity_VI=ulambda.*elambda_y + vteta.*eteta_y;
vitz_VI=ulambda.*elambda_z + vteta.*eteta_z;
end