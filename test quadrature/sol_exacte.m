function [vt] = sol_exacte(x,y,z)
global test
global radius

if test == 1
    [lambda, teta, ~]=cart2sph(x,y,z);
    xx=x./radius;
    yy=y./radius;
    zz=z./radius;
    
    uu=(xx).^2+xx.*yy.^3+6*cos(zz.^3);
    vv=sin(xx.*exp(-zz));

    elambda_x = -sin(lambda);
    elambda_y =  cos(lambda);
    elambda_z = zeros(size(x));
    eteta_x = -sin(teta).*cos(lambda);
    eteta_y = -sin(teta).*sin(lambda);
    eteta_z =  cos(teta);

    vt(:,:,1)=uu.*elambda_x + vv.*eteta_x;
    vt(:,:,2)=uu.*elambda_y + vv.*eteta_y;
    vt(:,:,3)=uu.*elambda_z + vv.*eteta_z;
    
end
    
end

