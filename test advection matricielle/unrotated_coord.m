function [ lambda, teta ] = unrotated_coord(lambda_prime, teta_prime)
% passage de coordonnées en rotations à cooronnées sphériques classiques
global lambda_p teta_p
nom=cos(teta_prime).*sin(lambda_prime);
denom=sin(teta_prime).*cos(teta_p) + cos(teta_prime).*cos(lambda_prime).*sin(teta_p);
lambda = lambda_p + atanp(nom,denom);
teta=asin(sin(teta_prime).*sin(teta_p)-cos(teta_prime).*cos(teta_p).*cos(lambda_prime));
end