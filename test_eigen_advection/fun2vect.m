function [ fun ] = fun2vect( funfI,funfII,funfIII,funfIV,funfV,funfVI )
global nn
fI=reshape(funfI(1:nn-1,2:nn-1),[],1);
fII=reshape(funfII(1:nn-1,2:nn-1),[],1);
fIII=reshape(funfIII(1:nn-1,2:nn-1),[],1);
fIV=reshape(funfIV(1:nn-1,2:nn-1),[],1);
fV=reshape(funfV,[],1);
fVI=reshape(funfVI,[],1);
fun=[fI; fII; fIII; fIV; fV; fVI];
end

