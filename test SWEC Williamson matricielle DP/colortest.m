figure(4)
hFig = figure(4);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 2000 500])
plot_cs102(n,nn,-vort_fI,-vort_fII,-vort_fIII,-vort_fIV,-vort_fV,-vort_fVI)
title(['vorticity at time : ', num2str(time(end))])
axis([-pi pi 0 pi/2])
colormap jet
caxis([-1.5*0 1]*10^-4)