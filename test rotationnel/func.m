function [ f, vort ] = func(x,y,z)
global radius
test=1;
if test == 1
    p=3;
    q=3;
    r=3;
    h=x.^p.*y.^q.*z.^r;
    [lambda, teta,~]=cart2sph(x,y,z);
    elx=-sin(lambda);
    ely=cos(lambda);
    elz=0;

    f(:,:,1)=h.*elx;
    f(:,:,2)=h.*ely;
    f(:,:,3)=h.*elz;

    dhteta=(p.*x.^(p-1).*y.^q.*z.^r).*(-radius.*sin(teta).*cos(lambda));
    dhteta=dhteta+(q.*x.^p.*y.^(q-1).*z.^r).*(-radius.*sin(teta).*sin(lambda));
    dhteta=dhteta+(r.*x.^p.*y.^q.*z.^(r-1)).*(radius.*cos(teta));


    vort=h.*tan(teta)./radius-(1./radius).*dhteta;
elseif test == 2
    p=1;
    q=2;
    r=3;
    h=x.^p.*y.^q.*z.^r;
    [lambda, teta,~]=cart2sph(x,y,z);
    etx=-sin(teta).*cos(lambda);
    ety=-sin(teta).*sin(lambda);
    etz=cos(teta);

    f(:,:,1)=h.*etx;
    f(:,:,2)=h.*ety;
    f(:,:,3)=h.*etz;
    
    dhx=p.*x.^(p-1).*y.^q.*z.^r;
    dhy=q.*x.^p.*y.^(q-1).*z.^r;
    
    vort=-dhx.*sin(lambda)+dhy.*cos(lambda);
    
end

