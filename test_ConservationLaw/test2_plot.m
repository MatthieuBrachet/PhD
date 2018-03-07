%% plot err test 2
if test == 2
    indI=(x_fI<=0);
    indII=(x_fII<=0);
    indIII=(x_fIII<=0);
    indIV=(x_fIV<=0);
    indV=(x_fV<=0);
    indVI=(x_fVI<=0);
else
    indI=ones(size(x_fI));
    indII=indI;
    indIII=indI;
    indIV=indI;
    indV=indI;
    indVI=indI;
end
err_fI=err_fI.*indI;
err_fII=err_fII.*indII;
err_fIII=err_fIII.*indIII;
err_fIV=err_fIV.*indIV;
err_fV=err_fV.*indV;
err_fVI=err_fVI.*indVI;

hFig = figure(100);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_cs102(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI);
colorbar
fig_placier