function [ lambda_prime, teta_prime ] = rotated_coord( lambda, teta )
% transformation des coordonnées sphériques dans le repère en rotation
global lambda_p teta_p
nom=cos(teta).*sin(lambda-lambda_p);
denom=cos(teta).*sin(teta_p).*cos(lambda-lambda_p)-cos(teta_p).*sin(teta);
lambda_prime = atanp(nom,denom);
teta_prime=asin(sin(teta).*sin(teta_p)+cos(teta).*cos(teta_p).*cos(lambda-lambda_p));
end