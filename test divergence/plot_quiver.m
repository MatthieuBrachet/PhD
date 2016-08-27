function []=plot_quiver(nn,vit_I,vit_II, vit_III, vit_IV, vit_V, vit_VI)
% PLOT VECTORS THE CUBED SPHERE GRID
% M. Brachet - SEPT 9, 2015

global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global radius;


%% PLOT OF A QUIVER ON THE SPHERE

% FACE F - I;
quiver3(x_fI, y_fI, z_fI,vit_I(:,:,1), vit_I(:,:,2), vit_I(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(radius,0,0,'F','Color','g','FontSize',15);

% FACE E - II;
quiver3(x_fII, y_fII, z_fII,vit_II(:,:,1), vit_II(:,:,2), vit_II(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,radius,0,'E','Color','g','FontSize',15);

% FACE B - III;
quiver3(x_fIII, y_fIII, z_fIII,vit_III(:,:,1), vit_III(:,:,2), vit_III(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(-radius,0,0,'B','Color','g','FontSize',15);

% FACE W - IV;
quiver3(x_fIV, y_fIV, z_fIV,vit_IV(:,:,1), vit_IV(:,:,2), vit_IV(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,-radius,0,'W','Color','g','FontSize',15);

% FACE N - V;
quiver3(x_fV, y_fV, z_fV,vit_V(:,:,1), vit_V(:,:,2), vit_V(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,0,radius,'N','Color','g','FontSize',15);

% FACE S - VI;
quiver3(x_fVI, y_fVI, z_fVI,vit_VI(:,:,1), vit_VI(:,:,2), vit_VI(:,:,3),'color',[1,0,0]);hold on;
axis([-radius radius -radius radius -radius radius]);
text(0,0,-radius,'S','Color','g','FontSize',15);

axis equal;
axis([-radius radius -radius radius -radius radius]);

end
