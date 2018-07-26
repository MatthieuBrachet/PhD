function [lambdai,coupe,cc] = coupe_lat(fun_fI,fun_fII,fun_fIII,fun_fIV,fun_fV,fun_fVI,lat)
global x_fI x_fII x_fIII x_fIV x_fV x_fVI
global y_fI y_fII y_fIII y_fIV y_fV y_fVI
global z_fI z_fII z_fIII z_fIV z_fV z_fVI

%% *** données
% **** points
[lambda_I,teta_I,~]=cart2sph(x_fI,y_fI,z_fI);
[lambda_II,teta_II,~]=cart2sph(x_fII,y_fII,z_fII);
[lambda_III,teta_III,~]=cart2sph(x_fIII,y_fIII,z_fIII);
[lambda_IV,teta_IV,~]=cart2sph(x_fIV,y_fIV,z_fIV);
[lambda_V,teta_V,~]=cart2sph(x_fV,y_fV,z_fV);
[lambda_VI,teta_VI,~]=cart2sph(x_fVI,y_fVI,z_fVI);

lambda=[lambda_I lambda_II lambda_III lambda_IV lambda_V lambda_VI];
teta=[teta_I teta_II teta_III teta_IV teta_V teta_VI];
fun=[fun_fI fun_fII fun_fIII fun_fIV fun_fV fun_fVI];

%% *** interpolation
x=linspace(-pi/2,pi/2,200);
y=lat.*ones(size(x));
[X,Y]=meshgrid(x,y);
Coupe = griddata(lambda,teta,fun,X,Y,'cubic');
coupe=Coupe(1,:);
lambdai=X(1,:);
cc=griddata(lambda,teta,fun,0,lat,'cubic');
end

