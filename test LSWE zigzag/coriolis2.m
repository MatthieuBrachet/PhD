function[flu_I,flu_II,flu_III,flu_IV,flu_V,flu_VI]=coriolis2(funfI,funfII,funfIII,funfIV,funfV,funfVI)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius omega



%% ------ face I
x=x_fI; y=y_fI; z=z_fI;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfI(:,:,1);
vel_y=funfI(:,:,2);
vel_z=funfI(:,:,3);
flu_I(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_I(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_I(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;

%% ------ face II
x=x_fII; y=y_fII; z=z_fII;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfII(:,:,1);
vel_y=funfII(:,:,2);
vel_z=funfII(:,:,3);
flu_II(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_II(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_II(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;

%% ------ face III
x=x_fIII; y=y_fIII; z=z_fIII;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfIII(:,:,1);
vel_y=funfIII(:,:,2);
vel_z=funfIII(:,:,3);
flu_III(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_III(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_III(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;

%% ------ face IV
x=x_fIV; y=y_fIV; z=z_fIV;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfIV(:,:,1);
vel_y=funfIV(:,:,2);
vel_z=funfIV(:,:,3);
flu_IV(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_IV(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_IV(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;

%% ------ face V
x=x_fV; y=y_fV; z=z_fV;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfV(:,:,1);
vel_y=funfV(:,:,2);
vel_z=funfV(:,:,3);
flu_V(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_V(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_V(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;

%% ------ face VI
x=x_fVI; y=y_fVI; z=z_fVI;
[~, teta,~]=cart2sph(x,y,z);
cste=2.*omega.*sin(teta)./radius;
fc_x(:,:)=cste.*x;
fc_y(:,:)=cste.*y;
fc_z(:,:)=cste.*z;
vel_x=funfVI(:,:,1);
vel_y=funfVI(:,:,2);
vel_z=funfVI(:,:,3);
flu_VI(:,:,1)= fc_y.*vel_z - fc_z.*vel_y;
flu_VI(:,:,2)= fc_z.*vel_x - fc_x.*vel_z;
flu_VI(:,:,3)= fc_x.*vel_y - fc_y.*vel_x;
