%%% Function f1 in Fornberg-Martel p.1176

function res = fun_f1( x,y,z )

    res = 1 + x + y.*y + x.*x.*y + x.*x.*x.*x +...
        y.*y.*y.*y.*y + x.*x.*y.*y.*z.*z;
    
end

