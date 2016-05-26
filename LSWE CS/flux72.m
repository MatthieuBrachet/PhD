function [flu_I,flu_II,flu_III,flu_IV,flu_V,flu_VI]=...
    flux72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global gp hp
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius omega



[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72(funfI(:,:,1:3),funfII(:,:,1:3),funfIII(:,:,1:3)...
        ,funfIV(:,:,1:3),funfV(:,:,1:3),funfVI(:,:,1:3),n,nn);
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr72(funfI(:,:,4),funfII(:,:,4),funfIII(:,:,4)...
        ,funfIV(:,:,4),funfV(:,:,4),funfVI(:,:,4),n,nn);

cste=2*omega/(radius*radius);

%% ------ face I
fc_x(:,:)=cste*x_fI.*z_fI;
fc_y(:,:)=cste*y_fI.*z_fI;
fc_z(:,:)=cste*z_fI.*z_fI;
vel_x=funfI(:,:,1);
vel_y=funfI(:,:,2);
vel_z=funfI(:,:,3);
flu_I(:,:,1)=gp*grad_I(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_I(:,:,2)=gp*grad_I(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_I(:,:,3)=gp*grad_I(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_I(:,:,4)=hp*div_fI;

%% ------ face II
fc_x(:,:)=cste*x_fII.*z_fII;
fc_y(:,:)=cste*y_fII.*z_fII;
fc_z(:,:)=cste*z_fII.*z_fII;
vel_x=funfII(:,:,1);
vel_y=funfII(:,:,2);
vel_z=funfII(:,:,3);
flu_II(:,:,1)=gp*grad_II(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_II(:,:,2)=gp*grad_II(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_II(:,:,3)=gp*grad_II(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_II(:,:,4)=hp*div_fII;

%% ------ face III
fc_x(:,:)=cste*x_fIII.*z_fIII;
fc_y(:,:)=cste*y_fIII.*z_fIII;
fc_z(:,:)=cste*z_fIII.*z_fIII;
vel_x=funfIII(:,:,1);
vel_y=funfIII(:,:,2);
vel_z=funfIII(:,:,3);
flu_III(:,:,1)=gp*grad_III(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_III(:,:,2)=gp*grad_III(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_III(:,:,3)=gp*grad_III(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_III(:,:,4)=hp*div_fIII;

%% ------ face IV
fc_x(:,:)=cste*x_fIV.*z_fIV;
fc_y(:,:)=cste*y_fIV.*z_fIV;
fc_z(:,:)=cste*z_fIV.*z_fIV;
vel_x=funfIV(:,:,1);
vel_y=funfIV(:,:,2);
vel_z=funfIV(:,:,3);
flu_IV(:,:,1)=gp*grad_IV(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_IV(:,:,2)=gp*grad_IV(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_IV(:,:,3)=gp*grad_IV(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_IV(:,:,4)=hp*div_fIV;

%% ------ face V
fc_x(:,:)=cste*x_fV.*z_fV;
fc_y(:,:)=cste*y_fV.*z_fV;
fc_z(:,:)=cste*z_fV.*z_fV;
vel_x=funfV(:,:,1);
vel_y=funfV(:,:,2);
vel_z=funfV(:,:,3);
flu_V(:,:,1)=gp*grad_V(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_V(:,:,2)=gp*grad_V(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_V(:,:,3)=gp*grad_V(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_V(:,:,4)=hp*div_fV;

%% ------ face VI
fc_x(:,:)=cste*x_fVI.*z_fVI;
fc_y(:,:)=cste*y_fVI.*z_fVI;
fc_z(:,:)=cste*z_fVI.*z_fVI;
vel_x=funfVI(:,:,1);
vel_y=funfVI(:,:,2);
vel_z=funfVI(:,:,3);
flu_VI(:,:,1)=gp*grad_VI(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
flu_VI(:,:,2)=gp*grad_VI(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
flu_VI(:,:,3)=gp*grad_VI(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
flu_VI(:,:,4)=hp*div_fVI;
