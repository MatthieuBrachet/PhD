% *************************************************************************
% test Linear Shallow Water
% 
% author :
%     - Matthieu Brachet
% *************************************************************************
clc; clear all; close all
%% *** global data ********************************************************
global n dxi radius u0;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test;

%% *** options ************************************************************
opt_ftr = 2;
test = 0; % (test=0 : personal test - 0 : divergence nulle; test=1 : galewsky regular; test=2 : galewsky par morceaux)
video = 'no';
ppp=4;

%% *** physical data ******************************************************

%% *** space data *********************************************************
n=40;
nn=n+2;
mod72;

%% *** time data **********************************************************
cfl=0.1;
u0=80;
ddt=cfl*radius*dxi/u0;
Tmax=12;
tmax=Tmax*86400;
itemax=0;
%% *** initialization *****************************************************
t=0; iter=0;
[ funfI] = fun (x_fI, y_fI, z_fI );
[ funfII ] = fun (x_fII, y_fII, z_fII );
[ funfIII ] = fun (x_fIII, y_fIII, z_fIII );
[ funfIV ] = fun (x_fIV, y_fIV, z_fIV );
[ funfV ] = fun (x_fV, y_fV, z_fV );
[ funfVI ] = fun (x_fVI, y_fVI, z_fVI );

%% *** video option *******************************************************
if strcmp(video,'yes')==1
    nFrames = min(itemax,floor(tmax/ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
end

%% time iteration
while t<tmax & iter<itemax
    clc;
    iter=iter+1
    t=t+ddt;
    
    %% filtrage
    for p=1:4
        [funfI(:,:,p),funfII(:,:,p),funfIII(:,:,p),funfIV(:,:,p),funfV(:,:,p),funfVI(:,:,p)]=...
            ftr72(funfI(:,:,p),funfII(:,:,p),funfIII(:,:,p),funfIV(:,:,p),funfV(:,:,p),funfVI(:,:,p),n,nn);
   end
    
    %% calcul de K0
    [kI,kII,kIII,kIV,kV,kVI]=flux72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
    
    %% Assemblage
    funfI=funfI-ddt*kI;
    funfII=funfII-ddt*kII;
    funfIII=funfIII-ddt*kIII;
    funfIV=funfIV-ddt*kIV;
    funfV=funfV-ddt*kV;
    funfVI=funfVI-ddt*kVI;
    
    %% film
    if strcmp(video,'yes')==1
        [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(funfI(:,:,1:3),funfII(:,:,1:3),funfIII(:,:,1:3)...
        ,funfIV(:,:,1:3),funfV(:,:,1:3),funfVI(:,:,1:3),n,nn);
        figure(2)
        plot_cs14(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI,-3*10^-4,10^-4)
        
%         figure(100)
%         plot_cs11(n,nn,funfI(:,:,ppp),funfII(:,:,ppp),funfIII(:,:,ppp),funfIV(:,:,ppp),funfV(:,:,ppp),funfVI(:,:,ppp))
        hold off;
        mov(iter) = getframe(gcf);
    end
end

ref=floor(10000*now);
if strcmp(video,'yes') == 1
    mkdir(['./video-' date ])
    movie2avi(mov, ['./video-' date '/ref_' num2str(ref) '.avi'], 'compression', 'None');
end

%% *** TRACÉ **************************************************************
figure(1)
plot_cs11(n,nn,funfI(:,:,ppp),funfII(:,:,ppp),funfIII(:,:,ppp),funfIV(:,:,ppp),funfV(:,:,ppp),funfVI(:,:,ppp))



[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(funfI(:,:,1:3),funfII(:,:,1:3),funfIII(:,:,1:3)...
        ,funfIV(:,:,1:3),funfV(:,:,1:3),funfVI(:,:,1:3),n,nn);
figure(2)
plot_cs11(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI)