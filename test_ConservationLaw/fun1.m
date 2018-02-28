function [fu] = fun1(x)
fu=1.*(x<=0)+(1-6*x.^2+(8/sqrt(2)).*x.^3).*(x>=0).*(x<=sqrt(2)/2);
end