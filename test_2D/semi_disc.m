function [kh,ku,kv] = semi_disc(h,u,v)
global gp hp coriolis
global px kx py ky
%% terme kh
w=kx*u;
dux=px\w;
w=ky*v;
dvy=py\w;
kh=-hp.*(dux+dvy);

%% terme ku
w=kx*h;
dhx=px\w;
ku=-gp.*dhx+coriolis.*v;

%% terme kv
w=ky*h;
dhy=py\w;
kv=-gp.*dhy-coriolis.*u;
end