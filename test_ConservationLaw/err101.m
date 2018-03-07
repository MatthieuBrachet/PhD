function [Er1,Er2,Eri] = err101(err_fI, err_fII, err_fIII,err_fIV,err_fV,err_fVI)
global n nn
global x_fI x_fII x_fIII
global x_fIV x_fV x_fVI
global test
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


str='1'; [~,~,~,~,~,~,Er1]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
str='2'; [~,~,~,~,~,~,Er2]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
str='infty'; [~,~,~,~,~,~,Eri]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
end

