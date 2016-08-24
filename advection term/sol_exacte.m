function [ht,vt] = sol_exacte(x,y,z,t)

global omega gp u0 h0 radius 
global teta0 teta1

    %% vitesse
    [n1,n2]=size(x);
    vt=zeros(n1,n2,3);
    [lambda, teta,r]=cart2sph(x,y,z);
    % uu=u0*fun10(teta,-pi/4, pi/4);
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
    %fa=radius*u0*fun10(a,-pi/4, pi/4).*(2*omega*sin(a));
    fa=radius*u0*fun10(a,teta0, teta1).*(2*omega*sin(a));
    %
    % fb=radius*u0*fun10(b,-pi/4, pi/4).*(2*omega*sin(b));
    fb=radius*u0*fun10(b,teta0, teta1).*(2*omega*sin(b));
    %
    int=0.5*(fa+fb).*h;
    for kk=1:n-1
        xx=a+kk*h;
        % fxx=radius*u0*fun10(xx,-pi/4, pi/4).*(2*omega*sin(xx));
        fxx=radius*u0*fun10(xx,teta0, teta1).*(2*omega*sin(xx));

        int=int+fxx.*h;
    end


    %% hauteur 4
    ht=h0-int/gp;
    
end