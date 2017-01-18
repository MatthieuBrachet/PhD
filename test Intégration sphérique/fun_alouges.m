function res = fun_alouges( X,x,y,z )

    ic=complex(0,1);
    res = exp(ic*(X(1)*x+X(2)*y+X(3)*z));

end

