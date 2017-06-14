function [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm_alpha(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,corner)

global dxi deta dga;


wei=ones(size(dga));
wei(1,1:end)=1/2;
wei(end,1:end)=1/2;
wei(1:end,1)=1/2;
wei(1:end,end)=1/2;
wei(1,1)=corner;
wei(1,end)=corner;
wei(end,1)=corner;
wei(end,end)=corner;
nrmI=dxi.*deta.*sum(sum(dga.*wei.*funfI));
nrmII=dxi.*deta.*sum(sum(dga.*wei.*funfII));
nrmIII=dxi.*deta.*sum(sum(dga.*wei.*funfIII));
nrmIV=dxi.*deta.*sum(sum(dga.*wei.*funfIV));
nrmV=dxi.*deta.*sum(sum(dga.*wei.*funfV));
nrmVI=dxi.*deta.*sum(sum(dga.*wei.*funfVI));
nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI;   
nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI;   
end