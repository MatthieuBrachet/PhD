function [ht,vt] = sol_exacte(x,y,z,t)
% exact solution for Williamson test.
global test
global gp u0 h0 radius omega
global alpha
global teta0 teta1

if test == -1
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;
    

elseif test == 0
    % test 2 of Williamson & al.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z=zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;

elseif test == 1
    % test 5 of Williamson & al.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;
    
elseif test == 2
    % test 5 of Williamson & al. wxith  smooth mountain.
    [lambda, teta, ~]=cart2sph(x,y,z);
    uu=u0.*(cos(teta).*cos(alpha)+cos(lambda).*sin(teta).*sin(alpha));
    vv=-u0.*sin(lambda).*sin(alpha);
    
    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;

    ht=h0-(1/gp)*(radius.*omega.*u0+0.5*u0.^2).*(-cos(lambda).*cos(teta).*sin(alpha)+sin(teta).*cos(alpha)).^2;
    
elseif test == 3
    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    uu=u0*fun10(teta,teta0, teta1);
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

    %% integrale

    nnn=10000;
    a=-pi/2;
    b=teta;
    fa=fun11(a);
    fb=fun11(b);
    h=(b-a)/nnn;

    s1=0;
    for kk=1:nnn/2-1
        x=a+(2*kk).*h;
        fx=fun11(x);
        s1=s1+fx;
    end

    s2=0;
    for kk=1:nnn/2
        x=a+(2*kk-1).*h;
        fx=fun11(x);
        s2=s2+fx;
    end

    int=(h/3).*(fa+2*s1+4*s2+fb);

    %% hauteur 4
    ht=h0-int/gp;
    
    
    
elseif test == 4
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
    
    nnn=10000;
    a=-pi/2;
    b=teta;
    fa=fun11(a);
    fb=fun11(b);
    h=(b-a)/nnn;

    s1=0;
    for kk=1:nnn/2-1
        x=a+(2*kk).*h;
        fx=fun11(x);
        s1=s1+fx;
    end

    s2=0;
    for kk=1:nnn/2
        x=a+(2*kk-1).*h;
        fx=fun11(x);
        s2=s2+fx;
    end

    int=(h/3).*(fa+2*s1+4*s2+fb);


    %% perturbation
    phihat=120;
    alphag=1/3;
    betag=1/15;
    teta2=pi/4;
    pert=phihat.*cos(teta).*exp(-(lambda./alphag).^2-((teta2-teta)./betag).^2);

    %% hauteur 4
    ht=h0-int/gp+pert;
end

