function c_mu = getRbfSphHarmCoef(mu,epsilon,rbf)
%   Returns the spherical harmonic expansion coefficients for the radial
%   kernel RBF. RBF is specified by a character string with the following
%   choices:
%
%      'MQ'  - multiquadric
%      'IMQ' - inverse multiquadric
%      'IQ'  - inverse quadratic
%      'GA'  - Gaussian

% Authors: Cecile Piret with minor modifications by Grady Wright
switch rbf
   case 'MQ'
      c_mu = -2*pi*(2*epsilon^2+1+(mu+1/2)*sqrt(1+4*epsilon^2))/...
      (mu+1/2)/(mu+3/2)/(mu-1/2)*(2/(1+sqrt(4*epsilon^2+1)))^(2*mu+1);
   case 'IMQ'
      c_mu = 4*pi/(mu+1/2)*(2/(1+sqrt(4*epsilon^2+1)))^(2*mu+1);
   case 'IQ'
      c_mu = 4*pi^(3/2)*factorial(mu)/gamma(mu+3/2)/(1+4*epsilon^2)^(mu+1)*...
      hypergeom([mu+1,mu+1],2*mu+2,4*epsilon^2/(1+4*epsilon^2));
   case 'GA'
      c_mu = 4*pi^(3/2)*exp(-2*epsilon^2)*besseli(mu+1/2,2*epsilon^2)/...
      epsilon^(2*mu+1);
end
