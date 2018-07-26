function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_xi101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[funftI1,funftII1,funftIII1,funftIV1,funftV1,funftVI1]=...
    ftr_eta101(funftI,funftII,funftIII,funftIV,funftV,funftVI,n,nn);

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_eta101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[funftI2,funftII2,funftIII2,funftIV2,funftV2,funftVI2]=...
    ftr_xi101(funftI,funftII,funftIII,funftIV,funftV,funftVI,n,nn);

funftI=.5*(funftI1+funftI2);
funftII=.5*(funftII1+funftII2);
funftIII=.5*(funftIII1+funftIII2);
funftIV=.5*(funftIV1+funftIV2);
funftV=.5*(funftV1+funftV2);
funftVI=.5*(funftVI1+funftVI2);
end

