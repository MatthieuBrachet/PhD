function [grade] = gr_harm(m,l,x,y,z)
[lambda, teta, r]=cart2sph(x,y,z);
phi=pi/2-teta;

[mm,pp]=size(x);
grade=zeros(mm,pp,3);

h_teta = dharmsph_theta(m,l,teta,lambda );
h_lambda = dharmsph_lambda(m,l,teta,lambda);

grade(:,:,1)=-1./r.*h_teta.*sin(teta) + 1./(r.*cos(teta)).*h_lambda.*cos(phi).*cos(teta);
grade(:,:,2)=1./r.*h_teta.*cos(teta) + 1./(r.*cos(teta)).*h_lambda.*cos(phi).*sin(teta);
grade(:,:,2)=1./(r.*sin(teta)).*h_lambda.*sin(phi);
end