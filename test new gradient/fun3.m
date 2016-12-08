function [v,v1,v2,v3]=fun3(x,y,z)
kw1=5;kw2=3;kw3=4;
v=sin(2*kw1*x*pi)+sin(2*kw2*y*pi)+sin(2*kw3*z*pi);
v1=2*kw1*pi*cos(2*kw1*x*pi);
v2=2*kw2*pi*cos(2*kw2*y*pi);
v3=2*kw3*pi*cos(2*kw3*z*pi);
end
