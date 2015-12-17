function s = rbfqrinterp(x,f,xi,rbf,epsilon)
%RBFQRINTERP Interpolates scattered data on the unit sphere using the RBF-QR 
% algorithm, which is suitable for "flat" RBFs (or small shape parameters).
%
%   S = rbfqrinterp(X,F,XI,RBF,EPSILON) Computes the RBF interpolant to the
%   the the target function F.  Returns the interpolant evaluated at the
%   points in XI.  Here X is assumed to be an N-by-3 array where each row
%   contains a point on the unit sphere in Cartesian coordiantes (x,y,z). N
%   should be a perfect square. F is assumed to be a N-by-k array where
%   each column contains samples of a different target function f at the
%   nodes in X.  XI is assumed to be an M-by-3 array where each row
%   contains a point on the unit sphere in Cartesian coordiantes where the
%   interpolant is to be evaluated. RBF is the radial kernel that is to be
%   used for constructing the interplant and is specified by a character
%   string with the following choices:
%     String  Name of kernel       Functional form of kernel phi(r)
%      'MQ'   multiquadric             sqrt(1 + (epsilon*r)^2)
%      'IMQ'  inverse multiquadric     1/sqrt(1 + (epsilon*r)^2)
%      'IQ'   inverse quadratic        1/(1 + (epsilon*r)^2)
%      'GA'   Gaussian                 exp(-(epsilon*r)^2)
%   Finally epsilon is the shape parameter that is to be used with the
%   radial kernel.
%   S will be an M-by-k vector where column j corresponds to the RBF
%   interpolant of the data in column j of F evaluated at the nodes in XI.
%
%   This function should only be used for 0 < Epsilon < 1.  For larger
%   shape parameters just use RBF-Direct.  Also, be aware that the function
%   can take considerable amount of time to complete, especially for large
%   N.
%
%   See the following paper for more details:
%      Fornberg B, Piret C (2007) A stable algorithm for flat radial basis
%      functions on a sphere. SIAM J. Sci. Comput. 30:60?80
%
%   Example:
%      f = @(x,y,z) x.*exp(y-z);  % Target function
%      xd = getMaxDetNodes(22^2); % Number nodes must be a perfect square
%      fd = f(xd(:,1),xd(:,2),xd(:,3));  % Target sampled at nodes x.
%      [X,Y,Z] = sphere(50);      % Where to evaluate the interpolant
%      xi = [X(:) Y(:) Z(:)];     % Arrange eval points as a M-by-3 array
%      rbf = 'IMQ';               % Use inverse multiquadric
%      epsilon = 0.05;            % Value of the shape parameter
%      s = rbfqrinterp(xd,fd,xi,rbf,epsilon);
%      s = reshape(s,size(X));    % Reshape interpolant for plotting 
%      surf(X,Y,Z,s);
%      shading interp; axis equal; colorbar;

% Authors: Cecile Piret with modifications by Grady Wright 2014

N = size(x,1);
if abs(round(sqrt(N))^2-N)/N > 10*eps
    error('RBFSPHERE:RBFQRINTERP','The number of nodes must be a perfect square');
end
% Switch to spherical coordinates.
[lambda,theta] = cart2sphm(x);

% Compute the new basis
[R,Y] = computeNewFlatBasis(lambda,theta,rbf,epsilon);

% Solve for the interpolation coefficients with respect to the new basis.
beta = Y*R\f;

% Evaluate the interpolant at the points specified.  This requires 
% switching the new basis.
[lambda_eval,theta_eval] = cart2sphm(xi);

% Each loop adds a block of columns of SPH of order mu to Y, evaluated at 
% the grid points
index = 1;
M = sqrt(size(R,1));
for mu = 0:M-1 
   Ye(:,index:2*mu+index) = blockSphHarm(mu,lambda_eval,theta_eval);
   index = index + 2*mu + 1;
end
% Compue the interpolant.
s = Ye*(R*beta);

end