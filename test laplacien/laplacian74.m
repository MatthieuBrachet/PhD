function [ A ] = laplacian74(n,nn)
for i=1:6*nn^2
    clc; disp([i 6*nn^2])
    X=zeros(6*nn^2,1);
    X(i)=1;
    [ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = deresh( X, n, nn );
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=gr74(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
    [lap_fI,lap_fII,lap_fIII,lap_fIV,lap_fV,lap_fVI]=div74(grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,n,nn);
    [ L ] = resh( lap_fI,lap_fII,lap_fIII,lap_fIV,lap_fV,lap_fVI, n, nn );
    A(:,i)=L;
end

A=sparse(A);