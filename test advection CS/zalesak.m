function [v] = zalesak(lambda, teta)
global radius
global lambdac1 tetac1
global lambdac2 tetac2
% calcul de la solution au temps t
RR=radius/3;
rdt1=radius*acos(sin(tetac1)*sin(teta)+cos(tetac1)*(cos(teta).*cos(lambda-lambdac1)));
rdt2=radius*acos(sin(tetac2)*sin(teta)+cos(tetac2)*(cos(teta).*cos(lambda-lambdac2)));
[n1,n2]=size(lambda);
% !!!! revoir le bloc 1 !!!!
for i=1:n1
    for j=1:n2
        if (rdt1(i,j) <= RR) & (abs(lambda(i,j) - lambdac1) >= RR/(6*radius))
            v(i,j) = 1;
        elseif ( rdt1(i,j) <= RR ) & ( abs(lambda(i,j)-lambdac1) <= RR/(6*radius) ) & (teta(i,j) - tetac1 < -5*RR/(12*radius))
            v(i,j) = 1;
        elseif (rdt2(i,j) <= RR) & (abs(lambda(i,j) - lambdac2) >= RR/(6*radius))
            v(i,j) = 1;
        elseif ( rdt2(i,j) <= RR ) & ( abs(lambda(i,j)-lambdac2) <= RR/(6*radius) ) & (teta(i,j) - tetac2 > 5*RR/(12*radius))
            v(i,j) = 1;
        else
            v(i,j) = 0.1;
        end
    end
end
end