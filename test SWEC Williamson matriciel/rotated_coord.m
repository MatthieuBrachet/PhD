function [ lambdap, tetap ] = rotated_coord( lambda, teta )
% transformation des coordonnées sphériques dans le repère en rotation
global alpha
tetap=asin(sin(teta).*cos(alpha)-cos(teta).*cos(lambda).*sin(alpha));
num=sin(teta)-sin(tetap).*cos(alpha);
denom=cos(tetap).*sin(alpha);
lp=acos(num./denom);
xwk=cos(alpha).*cos(lambda).*cos(teta)+sin(alpha).*sin(teta);
lambdap=lp.*(xwk<0)+(pi-lp).*(xwk>0);
