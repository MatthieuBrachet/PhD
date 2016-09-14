function [adv_fI,adv_fII,adv_fIII,adv_fIV,adv_fV,adv_fVI]=...
    adv73(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn)

%% *** terme 1 ************************************************************
% face I
norm2_I=vt_fI(:,:,1).^2+vt_fI(:,:,2).^2+vt_fI(:,:,3).^2;

% face II
norm2_II=vt_fII(:,:,1).^2+vt_fII(:,:,2).^2+vt_fII(:,:,3).^2;

% face III
norm2_III=vt_fIII(:,:,1).^2+vt_fIII(:,:,2).^2+vt_fIII(:,:,3).^2;

% face IV
norm2_IV=vt_fIV(:,:,1).^2+vt_fIV(:,:,2).^2+vt_fIV(:,:,3).^2;

% face V
norm2_V=vt_fV(:,:,1).^2+vt_fV(:,:,2).^2+vt_fV(:,:,3).^2;

% face VI
norm2_VI=vt_fVI(:,:,1).^2+vt_fVI(:,:,2).^2+vt_fVI(:,:,3).^2;

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr72(norm2_I, norm2_II, norm2_III, norm2_IV, norm2_V, norm2_VI,n,nn);

%% *** terme 2 ************************************************************

%% calcul du rotationnel par face
[rot_fI,rot_fII,rot_fIII,rot_fIV,rot_fV,rot_fVI]=...
    rot73(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

[n1,n2,n3]=size(rot_fI);
for i=1:n1
    for j=1:n2
        term_I(i,j,:)=cross(vt_fI(i,j,:),rot_fI(i,j,:));
        term_II(i,j,:)=cross(vt_fII(i,j,:),rot_fII(i,j,:));
        term_III(i,j,:)=cross(vt_fIII(i,j,:),rot_fIII(i,j,:));
        term_IV(i,j,:)=cross(vt_fIV(i,j,:),rot_fIV(i,j,:));
        term_V(i,j,:)=cross(vt_fV(i,j,:),rot_fV(i,j,:));
        term_VI(i,j,:)=cross(vt_fVI(i,j,:),rot_fVI(i,j,:));
    end
end

adv_fI   = 0.5 * grad_I   - term_I;
adv_fII  = 0.5 * grad_II  - term_II;
adv_fIII = 0.5 * grad_III - term_III;
adv_fIV  = 0.5 * grad_IV  - term_IV;
adv_fV   = 0.5 * grad_V   - term_V;
adv_fVI  = 0.5 * grad_VI  - term_VI;



end