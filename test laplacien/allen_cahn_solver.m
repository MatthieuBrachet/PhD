%% solveur Allen-Cahn on the CS.
clc; clear all; close all;

global n nn
global opt_ftr alfa_ftr scheme

video='no';
opt_ftr='redonnet10';
alfa_ftr=0;

%% *** numerical data *****************************************************
n=15; % for snapshot, n must be in the form 2^m-1 

scheme='explicite2';
mod74
[ MLAP ] = laplacian74(n,nn);
MLAP2=MLAP;

% scheme='compact4';
% mod74
% [ MLAP ] = laplacian74(n,nn);
% MLAP4=-MLAP;

ID=speye(size(MLAP2));
disp('matrix ok')

%% *** time data **********************************************************
Tmax=5;
ddt=0.01;
tau=1;
epsilon=0.1;

%% *** start problem ******************************************************
ht_fI=rand(nn,nn);
ht_fII=rand(nn,nn);
ht_fIII=rand(nn,nn);
ht_fIV=rand(nn,nn);
ht_fV=rand(nn,nn);
ht_fVI=rand(nn,nn);
U=resh(ht_fI,ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);
U=2*U-1;

%% *** time iterations ****************************************************
jour=date;
if strcmp(video,'yes')==1
    mkdir(['./RSS_video-' jour ])
    vidObj=VideoWriter(['./RSS_video-' jour 'AC.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

t=0;
while t<Tmax
    clc; t=t+ddt

    %% equations step
    % step 1 : solver for the heat equation
    Z=(ID+tau*ddt*MLAP2)\(-ddt*MLAP2*U);
    U=Z+U;
    
%     % step 2 : solver for non linear part
%     nom=U1;
%     denom=U1.^2+(1-U1.^2).*exp(-2*ddt./(epsilon.^2));
%     U=nom./sqrt(denom);
    
    
%     %% filtrage step
%     [ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = deresh( U, n, nn );
%     [funftI,funftII,funftIII,funftIV,funftV,funftVI]=ftr74(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
%     [ U ] = resh( funftI,funftII,funftIII,funftIV,funftV,funftVI, n, nn );
    
    %% video
    if strcmp(video,'yes')==1
        [ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = deresh( U, n, nn );
        
        figure(9)
        plot_cs7(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI);
        hold off
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
end

if strcmp(video,'yes') == 1
    close(vidObj);
end

%% *** figure *************************************************************
[ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = deresh( U, n, nn );

figure(1)
plot_cs11(n,nn,funfI,funfII, funfIII, funfIV, funfV, funfVI);


figure(2)
plot_cs7(n,nn,funfI,funfII, funfIII, funfIV, funfV, funfVI);

fig_placier