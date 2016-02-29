function [grade] = gr_harm(m,l,x,y,z)
[lambda, teta, r]=cart2sph(x,y,z);
phi=pi/2-teta;

[mm,pp]=size(x);
grade=zeros(mm,pp,3);

h_teta = dharmsph_theta(m,l,teta,lambda );
h_lambda = dharmsph_lambda(m,l,teta,lambda);

ii(:,:,1)=-sin(lambda);
ii(:,:,2)=cos(lambda);
ii(:,:,3)=zeros(size(lambda));

jj(:,:,1)=-sin(teta).*cos(lambda);
jj(:,:,2)=-sin(teta).*sin(lambda);
jj(:,:,3)=cos(teta);

grade(:,:,1)=1./(r.*cos(teta)).*h_lambda.*ii(:,:,1)+(1./r).*h_teta.*jj(:,:,1);
grade(:,:,2)=1./(r.*cos(teta)).*h_lambda.*ii(:,:,2)+(1./r).*h_teta.*jj(:,:,2);
grade(:,:,3)=1./(r.*cos(teta)).*h_lambda.*ii(:,:,3)+(1./r).*h_teta.*jj(:,:,3);
end