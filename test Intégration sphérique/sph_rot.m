function res = sph_rot( nhs,mhs,x,y,z,teta0,lambda0 )

    global xi;
    
    if nhs>0,
        [lambda,teta,radius]=cart2sph(x,y,z);
        ic=complex(0,1);
        nwk=legendre(nhs,sin(teta-teta0));
        mhs_a=abs(mhs);
        if size(x,2)==1
            nwk1=nwk(mhs_a+1,:)';
        else
            nwk1=squeeze(nwk(mhs_a+1,:,:));
        end
        Nlm=(-1)^mhs_a.*sqrt((2*nhs+1)/(4*pi)*factorial(nhs-mhs_a)/factorial(nhs+mhs_a));
        vc=...%1/sqrt(4*pi).*...
            Nlm.*nwk1.*exp(ic*mhs*(lambda-lambda0));
        res=vc;
    else
        res=1/sqrt(4*pi)*ones(size(x));
    end

end

