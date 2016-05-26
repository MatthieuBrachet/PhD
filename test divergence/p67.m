clear all;clc; close all;
global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
%%%%%%%%%%%%%%%%%%%

mod67; % APPEL MODULE "PROBLEME"

[mfunfI,dfunI]=fun6(x_fI,y_fI,z_fI); % FONCTION VECTORIELLE SUR LES 6 FACES
[mfunfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[mfunfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[mfunfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[mfunfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[mfunfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div68(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

% calcul erreur finale
for i=1:nn,
    for j=1:nn,
        errg_fI(i,j)=max(abs(div_fI(i,j)-dfunI(i,j)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fII(i,j)=max(abs(div_fII(i,j)-dfunII(i,j)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fIII(i,j)=max(abs(div_fIII(i,j)-dfunIII(i,j)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fIV(i,j)=max(abs(div_fIV(i,j,:)-dfunIV(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fV(i,j)=max(abs(div_fV(i,j,:)-dfunV(i,j,:)));
    end
end
for i=1:nn,
    for j=1:nn,
        errg_fVI(i,j)=max(abs(div_fVI(i,j,:)-dfunVI(i,j,:)));
    end
end
errI_25=max(max(abs(errg_fI)))
errII_26=max(max(abs(errg_fII)))
errIII_27=max(max(abs(errg_fIII)))
errIV_28=max(max(abs(errg_fIV)))
errV_29=max(max(abs(errg_fV)))
errVI_30=max(max(abs(errg_fVI)))
err_div_31=max([errVI_30,errV_29,errIV_28,errIII_27,errII_26,errI_25])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot_cs5(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI);colorbar;
title('DIV CALCULATED');
figure(2);
plot_cs5(n,nn,dfunI,dfunII,dfunIII,dfunIV,dfunV,dfunVI);colorbar;
title('DIV EXACT');
figure(3);
plot_cs5(n,nn,errg_fI,errg_fII,errg_fIII,errg_fIV,errg_fV,errg_fVI);colorbar;
title('ERREUR');

