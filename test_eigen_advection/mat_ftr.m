function [ sm ] = mat_ftr( vect )
global n nn
[ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = vect2fun( vect );
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[ sm ] = fun2vect(funftI,funftII,funftIII,funftIV,funftV,funftVI);
end