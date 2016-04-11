function [x] = SMW_tridiagper(a, b, f)
%**************************************************************************
%                     SOLVEUR TRIDIAGONAL PERIODIQUE
%
% Résolution de Ax=f
% avec A tridiagonale periodique symétrique a diagonale constante.
%
% avec A = [ a1   b1   0               0   bN]
%          [ b2   a2   b2   0              0 ]
%          [ 0    b3   a3   b3  0          0 ]
%          [           ...  ...  ...         ]
%          [ b1   0                    bN aN ]
%
% en utilisant la formule de Shermann Morisson Woodbury
% et un solveur tridiagonal de type algorithme de Thomas.
%
%**************************************************************************

n=length(a);

% matrices de periodicité
R=zeros(n,2);
R(1,1)=b(end);
R(end,end)=b(1);
S=zeros(n,2);
S(end,1)=1;
S(1,end)=1;

% matrice tridiagonale
B=sparse(diag(a)+diag(b(1,end-1),1)+diag(b(2:end),-1));
Id=speye(2,2);

% solveur de Shermann Morisson Woodbury
V1=B\f;
V2=(S')*V1;
mat=B\R;
mat=Id+S'*mat;
V3=mat\V2;
V4=R*V3;
V5=B\V4;
x=V1-V5;

end

