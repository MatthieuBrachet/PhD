function [ X ] = resh( funfI, funfII, funfIII, funfIV, funfV, funfVI, n, nn )
X1=reshape(funfI,[],1);
X2=reshape(funfII,[],1);
X3=reshape(funfIII,[],1);
X4=reshape(funfIV,[],1);
X5=reshape(funfV,[],1);
X6=reshape(funfVI,[],1);
X=[X1;X2;X3;X4;X5;X6];
end

