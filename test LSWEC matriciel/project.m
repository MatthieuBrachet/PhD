function [ v_lambda, v_teta ] = project(vt,x,y,z)
% Projection of a vector fields on the lambda component.
[lambda,teta,aaa]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
v_lambda=vt(:,:,1).*elambda_x+vt(:,:,2).*elambda_y+vt(:,:,3).*elambda_z;

eteta_x = -sin(teta).*cos(lambda);
eteta_y = -sin(teta).*sin(lambda);
eteta_z =  cos(teta);
v_teta=vt(:,:,1).*eteta_x+vt(:,:,2).*eteta_y+vt(:,:,3).*eteta_z;
end
