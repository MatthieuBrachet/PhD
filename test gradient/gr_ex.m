function [grade] = gr_ex(x,y,z)
[l,t,r]=cart2sph(x,y,z);
dudl=exp(l);
dudt=exp(t);

grade(:,:,1)=(1./(r.*cos(t))) .*dudl .*(-sin(l)) + (1./r).*dudt .*(-sin(t).*cos(l));
grade(:,:,2)=(1./(r.*cos(t))) .*dudl .*(cos(l))  + (1./r).*dudt .*(-sin(t).*sin(l));
grade(:,:,3)=(1./(r.*cos(t))) .*dudl .*(0)       + (1./r).*dudt .*(cos(t));
end