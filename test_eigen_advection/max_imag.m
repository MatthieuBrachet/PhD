function [ lambda ] = max_imag( V )
[~,i]=max(abs(imag(V)));
lambda=V(i);
end

