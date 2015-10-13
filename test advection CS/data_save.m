function []=data_save(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
% PLOT THE CUBED SPHERE GRID
% J-P. CROISILLE - JULY 20 2010
global test
global u0 radius
global lambda0 teta0

%% tracé de la courbe
r=radius;
%
xi=linspace(-pi/4, pi/4, nn); 
eta=linspace(-pi/4, pi/4, nn); 
for i=1:nn,
  for j=1:nn,
    xx(i,j)=tan(xi(i));
    yy(i,j)=tan(eta(j));
  end
end
%
for i=1:nn,
  for j=1:nn,
    delta(i,j)=1+xx(i,j)^2+yy(i,j)^2;
  end
end


fid = fopen('./results/data.csv','w');
fprintf(fid,'%s\n', ['x coord, y coord, z coord, scalar']);


% PLOT OF A POLYNOMIAL FUNCTION ON THE SPHERE
% FACE F - I;
x1=sqrt(radius^2./delta);
y1=xx.*x1;
z1=yy.*x1;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x1(i,j)) ', ', num2str(y1(i,j)) ', ', num2str(z1(i,j)) ', ', num2str(funfI(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');
        
     
% FACE E - II;
y2=sqrt(radius^2./delta);
x2=-xx.*y2;
z2=yy.*y2;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x2(i,j)) ', ', num2str(y2(i,j)) ', ', num2str(z2(i,j)) ', ', num2str(funfII(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');

% FACE B - III;
x3=-sqrt(radius^2./delta);
y3=xx.*x3;
z3=-yy.*x3;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x3(i,j)) ', ', num2str(y3(i,j)) ', ', num2str(z3(i,j)) ', ', num2str(funfIII(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');

% FACE W - IV;
y4=-sqrt(radius^2./delta);
x4=-xx.*y4;
z4=-yy.*y4;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x4(i,j)) ', ', num2str(y4(i,j)) ', ', num2str(z4(i,j)) ', ', num2str(funfIV(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');

% FACE N - V;
z5=sqrt(radius^2./delta);
y5=xx.*z5;
x5=-yy.*z5;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x5(i,j)) ', ', num2str(y5(i,j)) ', ', num2str(z5(i,j)) ', ', num2str(funfV(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');

% FACE S - VI;
z6=-sqrt(radius^2./delta);
y6=-xx.*z6;
x6=-yy.*z6;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        fprintf(fid,'%s\n',[num2str(x6(i,j)) ', ', num2str(y6(i,j)) ', ', num2str(z6(i,j)) ', ', num2str(funfVI(i,j))]);
    end
%    fprintf(fid,'%s\n','  ');
end
%fprintf(fid,'%s\n','  ');
%fprintf(fid,'%s\n','  ');

fclose(fid);

