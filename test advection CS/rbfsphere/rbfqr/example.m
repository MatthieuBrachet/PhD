% clear all; close all
epsilon = 0.9;
% Basis function; valid choices: 'MQ','IMQ','IQ','GA'
rbf = 'MQ'; 

% Nodes
N = (22)^2;
xd = getMaxDetNodes(N);

% Test function to interpolate
f = @(x,y,z) x.*exp(y-z); % Test function to interpolate
f = @(x,y,z) max(z,0);

% Evaluate the data values to interpolate
fd = f(xd(:,1),xd(:,2),xd(:,3));

% Resolution of the grid for evaluating the interpolant
res = 50;
[X,Y,Z] = sphere(res);
% Arrange these points for the rbfqrinterp function:
xi = [X(:) Y(:) Z(:)];

% Interpolate the data.
s = rbfqrinterp(xd,fd,xi,rbf,epsilon);

% Reshape the output to match the "meshgrid" arrangement of the evaluations
% ponts.  This is necessary for plotting with surf.
s = reshape(s,size(X));

%% Plotting
subplot(1,2,1); % plot the interpolant
surf(X,Y,Z,s);
shading interp;
axis equal;
colorbar

subplot(1,2,2); % plot the error
surf(X,Y,Z,s-f(X,Y,Z));
shading interp;
axis equal;
colorbar


