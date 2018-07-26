function [ lambda ] = max_real( V )
[~,i]=max(abs(real(V)));
lambda=V(i);
end

