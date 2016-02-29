function [mflu_I,mflu_II,mflu_III,mflu_IV,mflu_V,mflu_VI]=flu72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius gp hp omega
% FLUX EQUATIONS LTE
mflu_I=zeros(nn,nn,4);
mflu_II=zeros(nn,nn,4);
mflu_III=zeros(nn,nn,4);
mflu_IV=zeros(nn,nn,4);
mflu_V=zeros(nn,nn,4);
mflu_VI=zeros(nn,nn,4);
%
[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72(mfunfI(:,:,1:3),mfunfII(:,:,1:3),mfunfIII(:,:,1:3)...
    ,mfunfIV(:,:,1:3),mfunfV(:,:,1:3),mfunfVI(:,:,1:3),n,nn);
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr72(mfunfI(:,:,4),mfunfII(:,:,4),mfunfIII(:,:,4)...
    ,mfunfIV(:,:,4),mfunfV(:,:,4),mfunfVI(:,:,4),n,nn);
%
%gp=9.80616; % CSTE GRAVITATION
%hp=100; % MEAN HEIGTH
% omega= 6.37;  % ANGULAR ROTATION  %omega= 0.0d+00;  % ANGULAR ROTATION  
%
% norm_xI=x_fI/radius;norm_yI=y_fI/radius;norm_zI=z_fI/radius;
% norm_xII=x_fII/radius;norm_yII=y_fII/radius;norm_zII=z_fII/radius;
% norm_xIII=x_fIII/radius;norm_yIII=y_fIII/radius;norm_zIII=z_fIII/radius;
% norm_xIV=x_fIV/radius;norm_yIV=y_fIV/radius;norm_zIV=z_fIV/radius;
% norm_xV=x_fV/radius;norm_yV=y_fV/radius;norm_zV=z_fV/radius;
% norm_xVI=x_fVI/radius;norm_yVI=y_fVI/radius;norm_zVI=z_fVI/radius;
% %
cste=2*omega/(radius*radius);
% ------ face I
fc_x(:,:)=cste*x_fI.*z_fI;
fc_y(:,:)=cste*y_fI.*z_fI;
%fc_z(:,:)=cste*z_fI.*z_fI;
fc_z(:,:)=-cste*(x_fI.*x_fI+y_fI.*y_fI);
vel_x=mfunfI(:,:,1);
vel_y=mfunfI(:,:,2);
vel_z=mfunfI(:,:,3);
mflu_I(:,:,1)=gp*grad_I(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_I(:,:,2)=gp*grad_I(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_I(:,:,3)=gp*grad_I(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_I(:,:,4)=hp*div_fI;
% ------ face II
fc_x(:,:)=cste*x_fII.*z_fII;
fc_y(:,:)=cste*y_fII.*z_fII;
%fc_z(:,:)=cste*z_fII.*z_fII;
fc_z(:,:)=-cste*(x_fII.*x_fII+y_fII.*y_fII);
vel_x=mfunfII(:,:,1);
vel_y=mfunfII(:,:,2);
vel_z=mfunfII(:,:,3);
mflu_II(:,:,1)=gp*grad_II(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_II(:,:,2)=gp*grad_II(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_II(:,:,3)=gp*grad_II(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_II(:,:,4)=hp*div_fII;
% ------ face III
fc_x(:,:)=cste*x_fIII.*z_fIII;
fc_y(:,:)=cste*y_fIII.*z_fIII;
%fc_z(:,:)=cste*z_fIII.*z_fIII;
fc_z(:,:)=-cste*(x_fIII.*x_fIII+y_fIII.*y_fIII);
vel_x=mfunfIII(:,:,1);
vel_y=mfunfIII(:,:,2);
vel_z=mfunfIII(:,:,3);
mflu_III(:,:,1)=gp*grad_III(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_III(:,:,2)=gp*grad_III(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_III(:,:,3)=gp*grad_III(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_III(:,:,4)=hp*div_fIII;
% ------ face IV
fc_x(:,:)=cste*x_fIV.*z_fIV;
fc_y(:,:)=cste*y_fIV.*z_fIV;
%fc_z(:,:)=cste*z_fIV.*z_fIV;
fc_z(:,:)=-cste*(x_fIV.*x_fIV+y_fIV.*y_fIV);
vel_x=mfunfIV(:,:,1);
vel_y=mfunfIV(:,:,2);
vel_z=mfunfIV(:,:,3);
mflu_IV(:,:,1)=gp*grad_IV(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_IV(:,:,2)=gp*grad_IV(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_IV(:,:,3)=gp*grad_IV(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_IV(:,:,4)=hp*div_fIV;
% ------ face V
fc_x(:,:)=cste*x_fV.*z_fV;
fc_y(:,:)=cste*y_fV.*z_fV;
%fc_z(:,:)=cste*z_fV.*z_fV;
fc_z(:,:)=-cste*(x_fV.*x_fV+y_fV.*y_fV);
vel_x=mfunfV(:,:,1);
vel_y=mfunfV(:,:,2);
vel_z=mfunfV(:,:,3);
mflu_V(:,:,1)=gp*grad_V(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_V(:,:,2)=gp*grad_V(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_V(:,:,3)=gp*grad_V(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_V(:,:,4)=hp*div_fV;
% ------ face VI
fc_x(:,:)=cste*x_fVI.*z_fVI;
fc_y(:,:)=cste*y_fVI.*z_fVI;
%fc_z(:,:)=cste*z_fVI.*z_fVI;
fc_z(:,:)=-cste*(x_fVI.*x_fVI+y_fVI.*y_fVI);
vel_x=mfunfVI(:,:,1);
vel_y=mfunfVI(:,:,2);
vel_z=mfunfVI(:,:,3);
mflu_VI(:,:,1)=gp*grad_VI(:,:,1)+ fc_y.*vel_z - fc_z.*vel_y;
mflu_VI(:,:,2)=gp*grad_VI(:,:,2)+ fc_z.*vel_x - fc_x.*vel_z;
mflu_VI(:,:,3)=gp*grad_VI(:,:,3)+ fc_x.*vel_y - fc_y.*vel_x;
mflu_VI(:,:,4)=hp*div_fVI;
%--------