function [funftI1,funftII1,funftIII1,funftIV1,funftV1,funftVI1] = ftr_ns101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_xi101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[funftI1,funftII1,funftIII1,funftIV1,funftV1,funftVI1]=...
    ftr_eta101(funftI,funftII,funftIII,funftIV,funftV,funftVI,n,nn);
end

