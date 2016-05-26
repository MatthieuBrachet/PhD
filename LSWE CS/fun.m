function [ funf ] = fun ( x, y, z )
%% rmq pbm de singularité au pole lorsque nn est impaire!
global nn radius test

if test==0

    omega=7.292d-05
    hp=10000;
    gp=9.80616;

    %% vitesse 1:3
    n1=size(x,1);
    n2=size(x,2);
    funf=zeros(n1,n2,4);

    [lambda, teta,~]=cart2sph(x,y,z);
    
    uu=cos(teta);
    vv=zeros(size(uu));

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    funf(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    funf(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    funf(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    %% integrale
    n=200;
    a=0;
    b=teta;
    h=(b-a)/n;
    fa=radius*funfu(a).*(2*omega*sin(a));
    fb=radius*funfu(b).*(2*omega*sin(b));
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        fxx=radius*funfu(xx).*(2*omega*sin(xx));
        int=int+fxx.*h;
    end


    %% hauteur 4
    funf(:,:,4)=hp-int/gp;
    
    
else
    
    omega=7.292d-05;
    hp=10000;
    gp=9.80616;

    %% vitesse 1:3
    funf=zeros(nn,nn,4);
    [lambda, teta,~]=cart2sph(x,y,z);
    u=cos( teta );
    v=zeros(size(u));

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=0;
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    funf(:,:,1)=u.*elambda_x + v.*eteta_x;
    funf(:,:,2)=u.*elambda_y + v.*eteta_y;
    funf(:,:,3)=u.*elambda_z + v.*eteta_z;


    %% integrale
    n=200;
    a=0;
    b=teta;
    h=(b-a)/n;
    fa=radius*funfu(a).*(2*omega*sin(a)+tan(a)/radius.*funfu(a));
    fb=radius*funfu(b).*(2*omega*sin(b)+tan(b)/radius.*funfu(b));
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        fxx=radius*funfu(xx).*(2*omega*sin(xx)+tan(xx)/radius.*funfu(xx));
        int=int+fxx.*h;
    end

    %% perturbation
    perthat=120;
    alfa=3;
    beta=15;
    lambda0=pi/4;
    teta0=pi/4;
    pert=perthat.*sech(alfa.*(lambda-lambda0)).^2.*sech(beta.*(teta-teta0)).^2;

    %% hauteur 4
    funf(:,:,4)=hp-int/gp;%+pert;
end