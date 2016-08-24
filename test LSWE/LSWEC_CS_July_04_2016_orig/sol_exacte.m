function [ht,vt] = sol_exacte(x,y,z,t)

global test
global omega hp gp u0 h0 radius

if test == 0
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=u0*fun10(teta,-pi/4, pi/4);
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
    n=1000;
    a=0;
    b=teta;
    h=(b-a)/n;
    fa=radius*u0*fun10(a,-pi/4, pi/4).*(2*omega*sin(a));
    fb=radius*u0*fun10(b,-pi/4, pi/4).*(2*omega*sin(b));
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        fxx=radius*u0*fun10(xx,-pi/4, pi/4).*(2*omega*sin(xx));
        int=int+fxx.*h;
    end


    %% hauteur 4
    ht=h0-int/gp;
    
elseif test == 1
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=u0*cos( teta );
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

    %% hauteur
    cste=omega*radius*u0/gp;
    ht=h0-cste*sin(teta).^2;
    
elseif test == 2
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    pwr=5;
    uu=u0*cos( teta ).^pwr;
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

    %% hauteur
    ht=h0+2*radius*omega*u0/((pwr+1)*gp).*(cos(teta).^(pwr+1) -1);
   
elseif test == 3
    
    sigma=10^-4;
    
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    %teta=teta+pi;
    uu=(sqrt(gp*hp)/10).*cos(teta).^5.*exp(-sigma*t);
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
   ht=cos(teta).*sin(lambda).*exp(-sigma*t);

elseif test == 4
    sigma=10^-4;
    
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=(sqrt(gp*hp)/10).*cos(teta).*sin(sigma*t);
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
    pwa=1;
    pwb=5;
    ht=(cos(teta).^pwb).*(cos(3*lambda).^pwa).*sin(sigma*t);    
    
elseif test == 5
    sigma=10^-4;
    
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,~]=cart2sph(x,y,z);
    uu=(sqrt(gp*hp)/10).*fun10(teta,-pi/4,pi/4).*exp(-sigma*t);
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
    ht=fun10(teta,-pi/4,pi/4).*exp(-sigma*t);
  
    
end
end

