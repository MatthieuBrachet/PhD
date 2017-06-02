function [funftI,funftII,funftIII,funftIV,funftV,funftVI] = ftr_ns101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_xi101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_eta101(funftI,funftII,funftIII,funftIV,funftV,funftVI,n,nn);
end

