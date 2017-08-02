function [func,gr] = fun2(x,y,z)
global radius
test = 2;
if test == 1
    [lambda, teta, ~]=cart2sph(x,y,z);
    pwr=9;
    func=cos(teta).^pwr;

    elx=-sin(lambda);
    ely=cos(lambda);
    elz=0;

    etx=-sin(teta).*cos(lambda);
    ety=-sin(teta).*sin(lambda);
    etz=cos(teta);

    dhl=0;
    dht=-pwr*cos(teta).^(pwr-1).*sin(teta);

    gr(:,:,1)=1./(radius.*cos(teta)).*dhl.*elx+1./(radius).*dht.*etx;
    gr(:,:,2)=1./(radius.*cos(teta)).*dhl.*ely+1./(radius).*dht.*ety;
    gr(:,:,3)=1./(radius.*cos(teta)).*dhl.*elz+1./(radius).*dht.*etz;
    
    
    
elseif test == 2
    p=9;
    q=9;
    r=9;
    func=(x.^p).*(y.^q).*(z.^r);
    proj=(p+q+r)./radius.*func;
    
    gr(:,:,1)=p.*x.^(p-1).*y.^q.*z.^r-proj.*x./radius;
    gr(:,:,2)=q.*x.^p.*y.^(q-1).*z.^r-proj.*y./radius;
    gr(:,:,3)=r.*x.^p.*y.^q.*z.^(r-1)-proj.*z./radius;
end

end

