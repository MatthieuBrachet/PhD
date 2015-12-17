% ============================== FUNCTION RBFQR ====================================
function [R,Y] = computeNewFlatBasis(lambda,theta,rbf,epsilon)
%COMPUTENEWBASIS Computes a new basis for space spanned by shifts of "flat
% RBFs" on the sphere.
%
%   [R,Y] = computeNewFlatBasis(LAMBDA,THETA,RBF,EPSILON) Computes a new basis
%   for the space spanned by shifts of "flat RBFs" on the sphere, where the
%   shifts are specified by the spherical cooridnates (LAMBDA,THETA).  Here
%   LAMBDA is the azimuthal coordinate and THETA is the elevation
%   (meausured from the equator).  EPSILON is the value of the shape
%   parameter to use for the RBF.
%
%   RBF is specified by a character string with the following choices:
%
%      'MQ'  - multiquadric
%      'IMQ' - inverse multiquadric
%      'IQ'  - inverse quadratic
%      'GA'  - Gaussian
%
%   See the following paper for more details:
%      Fornberg B, Piret C (2007) A stable algorithm for flat radial basis
%      functions on a sphere. SIAM J. Sci. Comput. 30:60?80

% Authors: Cecile Piret with modifications by Grady Wright 2014

n = length(lambda); Y = zeros(n); B = zeros(n);
mu = 0; index = 1; orderDifference = 0;
mu_n = ceil(sqrt(n))-1; %the order of the n_th spherical harmonic
while orderDifference < -log10(eps) %eps is the machine precision
   % Each loop adds a block of columns of SPH of order mu to Y and to B.
   % Compute the spherical harmonics matrix
   Y(:,index:2*mu+index) = blockSphHarm(mu,lambda,theta);
   % Compute the expansion coefficients matrix
   B(:,index:2*mu+index) = Y(:,index:2*mu+index)*getRbfSphHarmCoef(mu,epsilon,rbf);
   B(:,index+mu) = B(:,index+mu)/2;
   % Truncation criterion
   if mu > mu_n-1
      orderDifference = log10(norm(B(:,mu_n^2+1:(mu_n+1)^2),inf)/...
      norm(B(:,(mu+1)^2),inf)*epsilon^(2*(mu_n-mu)));
   end
   index = index+2*mu+1; mu = mu+1; % Calculate column index of next block
end
[~,R] = qr(B); % QR-factorization to find the RBF_QR basis
E = epsilon.^(2*(repmat(ceil(sqrt(n+1:mu^2))-1,n,1) - ... % Introduce the
repmat(ceil(sqrt(1:n))-1,mu^2-n,1)')); % powers of epsilon
%Solve the interpolation linear system
R = [eye(n),E.*(R(1:n,1:n)\R(1:n,n+1:end))]';
