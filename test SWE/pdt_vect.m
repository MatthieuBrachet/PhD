function [w] = pdt_vect(x,y)
% pdt vectoriel de x par y (=cross(x,y)?)

w(1)=x(2)*y(3)-x(3)*y(2);
w(2)=x(3)*y(1)-x(1)*y(3);
w(3)=x(1)*y(2)-x(2)*y(1);
end

