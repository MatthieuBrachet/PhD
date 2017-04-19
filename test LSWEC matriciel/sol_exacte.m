function [ht,vt] = sol_exacte(x,y,z,t)
global n nn
global test
global omega hp gp u0 h0 radius
global teta0 teta1

if test == 0
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=u0*fun10(teta,teta0, teta1);
    vv=zeros(n1,n2);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    %% integrale
    nnn=1000;
    a=0;
    b=teta;
    h=(b-a)/nnn;
    fa=radius*u0*fun10(a,teta0, teta1).*(2*omega*sin(a));
    fb=radius*u0*fun10(b,teta0, teta1).*(2*omega*sin(b));
    int=0.5*(fa+fb).*h;
    for kk=1:nnn-1
        xx=a+kk*h;
        fxx=radius*u0*fun10(xx,teta0, teta1).*(2*omega*sin(xx));
        int=int+fxx.*h;
    end


    %% hauteur 4
    ht=h0-int/gp;

elseif test == 1
    sigma=10^-5;
    
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,aaa]=cart2sph(x,y,z);
    uu=(sqrt(gp*hp)/10).*fun10(teta,teta0,teta1).*exp(-sigma*t);
    vv=zeros(n1,n2);

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    %% hauteur
    ht=fun10(teta,teta0,teta1).*exp(-sigma*t);
elseif test == 2
    % Paldor test case (not ready).
    error('the Paldor test case is not ready yet.')
    [lambda, teta,aaa]=cart2sph(x,y,z);
    lev=5;
    waveFlag='EIG';
    [u,v,h]=shallow_water_waves_test(lambda,teta,lev,t,waveFlag);
    uu=u;
    vv=v;

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    %% hauteur
    ht=h;
    
elseif test == 3
    % personnal test case
    vt=zeros(nn,nn,3);
    [lambdat,tetat,~]=cart2sph(x,y,z);
    RR=radius/2;
    tetac=0;
    lambdac=0;
    rdt=radius*acos(sin(tetac)*sin(tetat)+cos(tetac)*(cos(tetat).*cos(lambdat-lambdac)));
    lwkt=(rdt<RR);
    for i=1:size(x,1)
        for j=1:size(x,2)  
            ht(i,j)=0.5*hp*(1+cos(pi*rdt(i,j)/RR))*lwkt(i,j);
        end
    end
end
end

