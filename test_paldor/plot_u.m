function [up_fI,up_fII,up_fIII,up_fIV,up_fV,up_fVI]=plot_u(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

%% face I
x=x_fI; y=y_fI; z=z_fI;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fI=vt_fI(:,:,1).*elambda_x+vt_fI(:,:,2).*elambda_y+vt_fI(:,:,3).*elambda_z;

%% face II
x=x_fII; y=y_fII; z=z_fII;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fII=vt_fII(:,:,1).*elambda_x+vt_fII(:,:,2).*elambda_y+vt_fII(:,:,3).*elambda_z;

%% face III
x=x_fIII; y=y_fIII; z=z_fIII;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fIII=vt_fIII(:,:,1).*elambda_x+vt_fIII(:,:,2).*elambda_y+vt_fIII(:,:,3).*elambda_z;

%% face IV
x=x_fIV; y=y_fIV; z=z_fIV;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fIV=vt_fIV(:,:,1).*elambda_x+vt_fIV(:,:,2).*elambda_y+vt_fIV(:,:,3).*elambda_z;

%% face V
x=x_fV; y=y_fV; z=z_fV;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fV=vt_fV(:,:,1).*elambda_x+vt_fV(:,:,2).*elambda_y+vt_fV(:,:,3).*elambda_z;

%% face VI
x=x_fVI; y=y_fVI; z=z_fVI;
[lambda,~,~]=cart2sph(x,y,z);
elambda_x = -sin(lambda);
elambda_y =  cos(lambda);
elambda_z = zeros(size(x));
up_fVI=vt_fVI(:,:,1).*elambda_x+vt_fVI(:,:,2).*elambda_y+vt_fVI(:,:,3).*elambda_z;
end

