function [ f, vort ] = func(x,y,z)
global radius
test=2;
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
    a=radius;
    p=4;
    q=2;
    r=3;
    h=(x/a).^p.*(y/a).^q.*(z/a).^r;
    [lambda, teta,~]=cart2sph(x,y,z);
    etx=-sin(teta).*cos(lambda);
    ety=-sin(teta).*sin(lambda);
    etz=cos(teta);

    f(:,:,1)=h.*etx;
    f(:,:,2)=h.*ety;
    f(:,:,3)=h.*etz;
    
    dhx=p.*(x/a).^(p-1).*(y/a).^q.*(z/a).^r;
    dhy=q.*(x/a).^p.*(y/a).^(q-1).*(z/a).^r;
    
    vort=-dhx.*sin(lambda)+dhy.*cos(lambda);
    vort=vort/a;
    
elseif test == 3
   f(:,:,1)=exp(y/radius);
   f(:,:,2)=exp(z/radius);
   f(:,:,3)=exp(x/radius);
   
   rx=-exp(z/radius)/radius;
   ry=-exp(x/radius)/radius;
   rz=-exp(y/radius)/radius;
   
   vort=(rx.*x+ry.*y+rz.*z)./radius;
   
elseif test == 4
    alpha=3;
    u0=1;
    [lambda, teta,~]=cart2sph(x,y,z);
    h=u0.*cos(teta).^alpha;
    elx=-sin(lambda);
    ely=cos(lambda);
    elz=0;

    f(:,:,1)=h.*elx;
    f(:,:,2)=h.*ely;
    f(:,:,3)=h.*elz;

    dh=-u0.*alpha.*sin(teta).*cos(teta).^(alpha-1);
    %vort=(1./radius).*(h.*tan(teta)-dh);
    
    vort=(alpha+1)/radius.*cos(teta).^(alpha-1).*sin(teta);
elseif test == 5
   a=radius;
   f(:,:,1)=y/a.*exp(y/a)-z/a.*exp(x/a);
   f(:,:,2)=z/a.*exp(z/a)-x/a.*exp(y/a);
   f(:,:,3)=x/a.*exp(x/a)-y/a.*exp(z/a);
   
   vort=-x/a.*(2+z/a).*exp(z/a)-y/a.*(2+x/a).*exp(x/a)-z/a.*(2+y/a).*exp(y/a);
   vort=vort./a;
end

