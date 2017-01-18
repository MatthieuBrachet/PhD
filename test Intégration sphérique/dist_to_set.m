function res = dist_to_set( point,Set )

    n=size(Set,1);
    res=Inf;
    for i=1:n,
        Xi=Set(i,:);
        tmp=real(acos(point*Xi'));
        if tmp<res,
            res=tmp;
        end
    end    

end

