function SPHBlockMu = blockSphHarm(mu,th,lam)
%SPHBLOCKMU Returns a matrix containing all spherical harmonics of degree 
% mu evaluated at the locations (lam,th).
N = length(th);
L_mu_nu(:,1:mu+1) = legendre(mu,sin(lam))'; 
a = 0:mu;
scale = repmat(sqrt(factorial(1+mu-a-1)./factorial(1+mu+a-1)),N,1);
Y = scale.*L_mu_nu(:,a+1).*exp(1i*repmat(a,N,1).*repmat(th,1,mu+1));
SPHBlockMu = sqrt((2*mu+1)/(4*pi))*[imag(Y(:,end:-1:2)),real(Y)];
end
