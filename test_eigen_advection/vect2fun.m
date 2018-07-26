function [ funI, funII, funIII, funIV, funV, funVI ] = vect2fun( fun )
global nn

Nb=(nn-2)*(nn-1);
fI=reshape(fun(1:Nb),nn-1,nn-2);
fII=reshape(fun(Nb+1:2*Nb),nn-1,nn-2);
fIII=reshape(fun(2*Nb+1:3*Nb),nn-1,nn-2);
fIV=reshape(fun(3*Nb+1:4*Nb),nn-1,nn-2);
ffV=reshape(fun(4*Nb+1:4*Nb+nn*nn),nn,nn);
ffVI=reshape(fun(4*Nb+nn*nn+1:end),nn,nn);

ffI=[fI; fII(1,:)];
ffII=[fII; fIII(1,:)];
ffIII=[fIII; fIV(1,:)];
ffIV=[fIV; fI(1,:)];

ffI=[ffVI(:,nn) ffI ffV(:,1)];
ffII=[ffVI(nn,nn:-1:1)' ffII ffV(nn,:)'];
ffIII=[ffVI(nn:-1:1,1) ffIII ffV(nn:-1:1,nn)];
ffIV=[ffVI(1,1:nn)' ffIV ffV(1,nn:-1:1)'];

funI=ffI;
funII=ffII;
funIII=ffIII;
funIV=ffIV;
funV=ffV;
funVI=ffVI;
end

