function [ y ] = tanc( x )
seuil=0.1;
y=tan(x).*(abs(x-pi/2)>seuil).*(abs(x-pi/2)>seuil);
end

