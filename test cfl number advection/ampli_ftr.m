function [ ampli ] = ampli_ftr( J, teta )
A=zeros(J+1,J+1);
NN=[0:J];
A(1,:)=1;
A(2,:)=(-1).^NN;
for i=3:J+1
    for j=2:J+1
        A(i,j)=(j-1).^(2*(i-2));
    end
end
b=zeros(J+1,1);
b(1)=1;
P=tril(A);
Ap=P\A;
bp=P\b;
a=Ap\bp;
%[a,FLAG,RELRES,ITER]=gmres(A,b);

ampli=zeros(size(teta));
for k=1:length(a)
    ampli=ampli+a(k).*cos((k-1)*teta);
end
ampli=abs(ampli);
end

