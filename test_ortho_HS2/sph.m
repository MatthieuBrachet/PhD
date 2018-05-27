function res = sph( nhs,mhs,x,y,z )
% function res = sph( nhs,mhs,x,y,z ).
% Spherical harmonic of order nhs,mhs ( |mhs| <= nhs )  
if (abs(mhs)>nhs)
    error('We must have  |mhs| <= nhs.')
end
if nhs>0,
    [lambda,teta,~]=cart2sph(x,y,z);
    ic=complex(0,1);
    nwk=legendre(nhs,sin(teta));
    mhs_a=abs(mhs);
    if size(x,2)==1
        nwk1=nwk(mhs_a+1,:)';
    else
        nwk1=squeeze(nwk(mhs_a+1,:,:));
    end
    Nlm=(-1)^mhs_a.*sqrt((2*nhs+1)/(4*pi)*factorial(nhs-mhs_a)/factorial(nhs+mhs_a));
    vc=Nlm.*nwk1.*exp(ic*mhs*lambda);
    res=vc;
else
    res=1/sqrt(4*pi)*ones(size(x));
end
end