function [ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = deresh( X, n, nn )
N=nn*nn;
funfI=reshape(X(1:N),nn,nn);
funfII=reshape(X(N+1:2*N),nn,nn);
funfIII=reshape(X(2*N+1:3*N),nn,nn);
funfIV=reshape(X(3*N+1:4*N),nn,nn);
funfV=reshape(X(4*N+1:5*N),nn,nn);
funfVI=reshape(X(5*N+1:6*N),nn,nn);
end