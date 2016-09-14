function [ht,vt] = sol_exacte(x,y,z,t)

global test
global gp u0 h0
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
    n=1000;
    a=0;
    b=teta;
    h=(b-a)/n;
    %fa
    fa=fun11(a);
    %
    % fb
    fb=fun11(b);
    %
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        fxx=fun11(xx);
        int=int+fxx.*h;
    end


    %% hauteur 4
    ht=h0-int/gp;
    
elseif test == 1
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
    n=1000;
    a=0;
    b=teta;
    h=(b-a)/n;
    %fa
    fa=fun11(a);
    %
    % fb
    fb=fun11(b);
    %
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        fxx=fun11(xx);
        int=int+fxx.*h;
    end

    %% perturbation
    phihat=120;
    alpha=1/3;
    beta=1/15;
    teta2=pi/4;
    pert=phihat.*cos(teta).*exp(-(lambda./alpha).^2).*exp(-((teta2-teta)./beta).^2);

    %% hauteur 4
    ht=h0-int/gp+pert;
    
end

