function [y] = fun(x)
% xx=x*2*pi;
% y=cos(xx).*sin(2*xx)+sin(xx);
% y=y./sqrt(2);
xx=x-floor(x);
y=(abs(xx-.5)<.25)-.5;
y=2*y;

end

